/*
    Copyright 2005-2007 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#include "tbb/concurrent_vector.h"
#include "tbb_misc.h"
#include <stdexcept>
#include "itt_notify.h"
#include "tbb/task.h"


namespace tbb {

namespace internal {
class concurrent_vector_base::helper {
public:
    inline static size_type find_segment_end(const concurrent_vector_base &v) {
        segment_index_t u = v.my_segment==(&(v.my_storage[0])) ? pointers_per_short_table
                                                               : pointers_per_long_table;
        segment_index_t k = 0;
        while( k < u && v.my_segment[k].array )
            ++k;
        return k;
    }

    static void assign_first_segment_if_neccessary(concurrent_vector_base &v, segment_index_t k) {
        if( !v.my_first_block ) {
            v.my_first_block.compare_and_swap(k+1, 0);
        }
    }

    inline static void *allocate_segment(concurrent_vector_base &v, size_type n) {
        return v.vector_allocator_ptr(v, n);
    }

    inline static size_type enable_segment(concurrent_vector_base &v, size_type k, size_type element_size) {
        __TBB_ASSERT( !v.my_segment[k].array, "concurrent operation during growth?" );
        size_type m = segment_size(k);
        if( !k ) {
            assign_first_segment_if_neccessary(v, default_initial_segments-1);
            v.my_segment[0].array = allocate_segment(v, segment_size(v.my_first_block) );
            return 2;
        }
        if( !v.my_first_block )
            internal::SpinwaitWhileEq( v.my_first_block, segment_index_t(0) );
        if( k < v.my_first_block ) {
            if( !v.my_segment[0].array )
                internal::SpinwaitWhileEq( v.my_segment[0].array, (void*)0 );
            v.my_segment[k].array = reinterpret_cast<void*>(
                reinterpret_cast<char*>(v.my_segment[0].array) + segment_base(k)*element_size );
        } else
            v.my_segment[k].array = allocate_segment(v, m);
        return m;
    }

    inline static void extend_table_if_necessary(concurrent_vector_base &v, size_type k) {
        if(k >= pointers_per_short_table && v.my_segment == v.my_storage)
            extend_segment_table(v);
    }

    static void extend_segment_table(concurrent_vector_base &v) {
        segment_t* s = (segment_t*)NFS_Allocate( pointers_per_long_table, sizeof(segment_t), NULL );
        memset( s, 0, pointers_per_long_table*sizeof(segment_t) );
        // If other threads are trying to set pointers in the short segment, wait for them to finish their
        // assigments before we copy the short segment to the long segment.
        ExponentialBackoff backoff;
        segment_index_t i = 0; do {
            if(!v.my_storage[i].array) {
                backoff.pause(); continue;
            } else i++;
        } while( i < pointers_per_short_table);

        for( segment_index_t i = 0; i < pointers_per_short_table; i++)
            s[i] = v.my_storage[i];
        if( v.my_segment.compare_and_swap( s, v.my_storage ) != v.my_storage )
            NFS_Free( s );
    }
};

concurrent_vector_base::~concurrent_vector_base() {
    segment_t* s = my_segment;
    if( s != my_storage ) {
        // Clear short segment.
        for( segment_index_t i = 0; i < pointers_per_short_table; i++)
            my_storage[i].array = NULL;
        my_segment = my_storage;
        NFS_Free( s );
    }
}

concurrent_vector_base::size_type concurrent_vector_base::internal_capacity() const {
    return segment_base( helper::find_segment_end(*this) );
}

void concurrent_vector_base::internal_throw_exception(size_type) const {
    throw std::out_of_range("Index out of range");
}

void concurrent_vector_base::internal_reserve( size_type n, size_type element_size, size_type max_size ) {
    if( n>max_size ) {
        throw std::length_error("argument to ConcurrentVector::reserve exceeds ConcurrentVector::max_size()");
    }
    helper::assign_first_segment_if_neccessary(*this, segment_index_of(n));
    for( segment_index_t k = helper::find_segment_end(*this); segment_base(k)<n; ++k ) {
        helper::extend_table_if_necessary(*this, k);
        helper::enable_segment(*this, k, element_size);
    }
}

void concurrent_vector_base::internal_copy( const concurrent_vector_base& src, size_type element_size, internal_array_op2 copy ) {
    size_type n = src.my_early_size;
    my_early_size = n;
    my_segment = my_storage;
    if( n ) {
        helper::assign_first_segment_if_neccessary(*this, segment_index_of(n));
        size_type b;
        for( segment_index_t k=0; (b=segment_base(k))<n; ++k ) {
            helper::extend_table_if_necessary(*this, k);
            size_type m = helper::enable_segment(*this, k, element_size);
            if( m > n-b ) m = n-b; 
            copy( my_segment[k].array, src.my_segment[k].array, m );
        }
    }
}

void concurrent_vector_base::internal_assign( const concurrent_vector_base& src, size_type element_size, internal_array_op1 destroy, internal_array_op2 assign, internal_array_op2 copy ) {
    size_type n = src.my_early_size;
    while( my_early_size>n ) { 
        segment_index_t k = segment_index_of( my_early_size-1 );
        size_type b=segment_base(k);
        size_type new_end = b>=n ? b : n;
        __TBB_ASSERT( my_early_size>new_end, NULL );
        destroy( (char*)my_segment[k].array+element_size*(new_end-b), my_early_size-new_end );
        my_early_size = new_end;
    }
    size_type dst_initialized_size = my_early_size;
    my_early_size = n;
    helper::assign_first_segment_if_neccessary(*this, segment_index_of(n));
    size_type b;
    for( segment_index_t k=0; (b=segment_base(k))<n; ++k ) {
        helper::extend_table_if_necessary(*this, k);
        if(!my_segment[k].array)
            helper::enable_segment(*this, k, element_size);
        size_type m = k? segment_size(k) : 2;
        if( m > n-b ) m = n-b;
        size_type a = 0;
        if( dst_initialized_size>b ) {
            a = dst_initialized_size-b;
            if( a>m ) a = m;
            assign( my_segment[k].array, src.my_segment[k].array, a );
            m -= a;
            a *= element_size;
        }
        if( m>0 )
            copy( (char*)my_segment[k].array+a, (char*)src.my_segment[k].array+a, m );
    }
    __TBB_ASSERT( src.my_early_size==n, "detected use of ConcurrentVector::operator= with right side that was concurrently modified" );
}

void* concurrent_vector_base::internal_push_back( size_type element_size, size_type& index ) {
    __TBB_ASSERT( sizeof(my_early_size)==sizeof(reference_count), NULL );
    size_type tmp = __TBB_FetchAndIncrementWacquire((tbb::internal::reference_count*)&my_early_size);
    index = tmp;
    segment_index_t k_old = segment_index_of( tmp );
    size_type base = segment_base(k_old);
    helper::extend_table_if_necessary(*this, k_old);
    segment_t& s = my_segment[k_old];
    if( !s.array ) {
        if( base==tmp ) {
            helper::enable_segment(*this, k_old, element_size);
            ITT_NOTIFY( sync_releasing, &s.array );
        } else {
            ITT_NOTIFY(sync_prepare, &s.array);
            internal::SpinwaitWhileEq( s.array, (void*)0 );
            ITT_NOTIFY(sync_acquired, &s.array);
        }
    }
    size_type j_begin = tmp-base;
    return (void*)((char*)s.array+element_size*j_begin);
}

void concurrent_vector_base::internal_grow_to_at_least( size_type new_size, size_type element_size, internal_array_op2 init, const void *src ) {
    size_type e = my_early_size;
    while( e<new_size ) {
        size_type f = my_early_size.compare_and_swap(new_size,e);
        if( f==e ) {
            internal_grow( e, new_size, element_size, init, src );
            return;
        }
        e = f;
    }
}

concurrent_vector_base::size_type concurrent_vector_base::internal_grow_by( size_type delta, size_type element_size, internal_array_op2 init, const void *src ) {
    size_type result = my_early_size.fetch_and_add(delta);
    internal_grow( result, result+delta, element_size, init, src );
    return result;
}

void concurrent_vector_base::internal_grow( const size_type start, size_type finish, size_type element_size, internal_array_op2 init, const void *src ) {
    __TBB_ASSERT( start<finish, "start must be less than finish" );
    size_type tmp = start;
    helper::assign_first_segment_if_neccessary(*this, segment_index_of(finish));
    do {
        segment_index_t k_old = segment_index_of( tmp );
        size_type base = segment_base(k_old);
        helper::extend_table_if_necessary(*this, k_old);
        segment_t& s = my_segment[k_old];
        if( !s.array ) {
            if( base==tmp ) {
                helper::enable_segment(*this, k_old, element_size);
                ITT_NOTIFY( sync_releasing, &s.array );
            } else {
                ITT_NOTIFY(sync_prepare, &s.array);
                internal::SpinwaitWhileEq( s.array, (void*)0 );
                ITT_NOTIFY(sync_acquired, &s.array);
            }
        }
        size_type n = k_old?segment_size(k_old):2;
        size_type j_begin = tmp-base;
        size_type j_end = n > finish-base ? finish-base : n;
        init( (void*)((char*)s.array+element_size*j_begin), src, j_end-j_begin );
        tmp = base+j_end;
    } while( tmp<finish );
}

concurrent_vector_base::segment_index_t concurrent_vector_base::internal_clear( internal_array_op1 destroy ) {
    // Set "my_early_size" early, so that subscripting errors can be caught.
    // FIXME - doing so may be hurting exception saftey
    __TBB_ASSERT( my_segment, NULL );
    size_type finish = my_early_size;
    my_early_size = 0;
    while( finish>0 ) {
        segment_index_t k_old = segment_index_of(finish-1);
        segment_t& s = my_segment[k_old];
        __TBB_ASSERT( s.array, NULL );
        size_type base = segment_base(k_old);
        size_type j_end = finish-base;
        __TBB_ASSERT( j_end, NULL );
        destroy( s.array, j_end );
        finish = base;
    }
    return helper::find_segment_end(*this);
}

concurrent_vector_base::segment_t *concurrent_vector_base::internal_compact( size_type element_size, void *table_space, internal_array_op1 destroy, internal_array_op2 copy )
{
    // TODO: #ifdef TBB_VECTOR_SECTIONS, free garbage
    static const size_type page_size = 4096;
    const segment_index_t k_stop = helper::find_segment_end(*this);
    if( !k_stop ) return NULL;
    const segment_index_t first_block = my_first_block;
    segment_index_t k = first_block;
    // TODO: consider on op: '%' or '<'
    while (k < k_stop && segment_size( k ) * element_size % page_size) k++;
    if ( k == first_block )
        return NULL;
    // start optimization
    my_first_block = k;
    memcpy(table_space, my_segment, k * sizeof(segment_t));
    segment_t *table = reinterpret_cast<segment_t*>(table_space);
    void *seg = helper::allocate_segment( *this, segment_size(k) );
    for (segment_index_t i = 0; i < k; i++) {
        void *s = my_segment[i].array = reinterpret_cast<void*>(
            reinterpret_cast<char*>(seg) + segment_base(i)*element_size );
        if( !i || i >= first_block) {
            size_type my_segment_size;
            if (!i) my_segment_size = segment_size( first_block );
            else {
                my_segment_size = segment_size( i );
                if (my_segment_size*2 > my_early_size)
                    my_segment_size -= my_segment_size*2 - my_early_size;
            }
            copy( s, table[i].array, my_segment_size );
            destroy( table[i].array, my_segment_size );
        }
    }
    return table;
}

void concurrent_vector_base::internal_swap(concurrent_vector_base& v)
{
    size_type my_sz = my_early_size, v_sz = v.my_early_size;
    if(!my_sz && !v_sz) return;
    my_early_size = v_sz; v.my_early_size = my_sz;
    size_type tmp = my_first_block; my_first_block = v.my_first_block; v.my_first_block = tmp;
    bool my_short = (my_segment == my_storage), v_short  = (v.my_segment == v.my_storage);
    if ( my_short && v_short ) { // swap both tables
        segment_t tbl[pointers_per_short_table];
        memcpy(tbl, my_storage, pointers_per_short_table * sizeof(segment_t));
        memcpy(my_storage, v.my_storage, pointers_per_short_table * sizeof(segment_t));
        memcpy(v.my_storage, tbl, pointers_per_short_table * sizeof(segment_t));
    }
    else if ( my_short ) { // my -> v
        memcpy(v.my_storage, my_storage, pointers_per_short_table * sizeof(segment_t));
        my_segment = v.my_segment; v.my_segment = v.my_storage;
    }
    else if ( v_short ) { // v -> my
        memcpy(my_storage, v.my_storage, pointers_per_short_table * sizeof(segment_t));
        v.my_segment = my_segment; my_segment = my_storage;
    } else {
        segment_t *ptr = my_segment; my_segment = v.my_segment; v.my_segment = ptr;
    }
}


} // namespace internal

} // tbb
