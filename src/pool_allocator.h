/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006
// Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.

/*
 * Copyright (c) 1996-1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */

#ifndef PIRANHA_POOL_ALLOCATOR_H
#define PIRANHA_POOL_ALLOCATOR_H

#include <bits/c++config.h>
#include <boost/static_assert.hpp>
#include <cstdlib>
#include <new>
#include <bits/functexcept.h>
#include <ext/atomicity.h>
#include <ext/concurrence.h>
#include <mm_malloc.h>

namespace
{
  __gnu_cxx::__mutex palloc_init_mutex;
}

namespace piranha
{

  using std::size_t;
  using std::ptrdiff_t;
  using __gnu_cxx::__atomic_add_dispatch;
  using __gnu_cxx::__mutex;
  using __gnu_cxx::__scoped_lock;

/**
 *  @brief  Base class for pool_allocator.
 *
 *  @if maint
 *  Uses various allocators to fulfill underlying requests (and makes as
 *  few requests as possible when in default high-speed pool mode).
 *
 *  Important implementation properties:
 *  0. If globally mandated, then allocate objects from new
 *  1. If the clients request an object of size > _S_max_bytes, the resulting
 *     object will be obtained directly from new
 *  2. In all other cases, we allocate an object of size exactly
 *     _S_round_up(requested_size).  Thus the client has enough size
 *     information that we can return the object to the proper free list
 *     without permanently losing part of the object.
 *
 *  @endif
 */
  template <int Alignment>
    class pool_allocator_base
  {
// Check that Alignment is > 0 and a multiple of 8
      BOOST_STATIC_ASSERT(Alignment > 0 && (Alignment%8) == 0);
    protected:

      enum { _S_align = Alignment };
      enum { _S_max_bytes = 128 };
      enum { _S_free_list_size = (size_t)_S_max_bytes / (size_t)_S_align };

      union _Obj
      {
        union _Obj* _M_free_list_link;
        char        _M_client_data[1];            // The client sees this.
      };

      static _Obj* volatile         _S_free_list[_S_free_list_size];

// Chunk allocation state.
      static char*                  _S_start_free;
      static char*                  _S_end_free;
      static size_t                 _S_heap_size;

      size_t
        _M_round_up(size_t __bytes)
        { return ((__bytes + (size_t)_S_align - 1) & ~((size_t)_S_align - 1)); }

      _Obj* volatile*
        _M_get_free_list(size_t __bytes);

      __mutex&
        _M_get_mutex();

// Returns an object of size __n, and optionally adds to size __n
// free list.
      void*
        _M_refill(size_t __n);

// Allocates a chunk for nobjs of size size.  nobjs may be reduced
// if it is inconvenient to allocate the requested number.
      char*
        _M_allocate_chunk(size_t __n, int& __nobjs);
  };

// Definitions for pool_allocator_base.
  template <int Alignment>
    inline typename pool_allocator_base<Alignment>::_Obj* volatile*
    pool_allocator_base<Alignment>::_M_get_free_list(size_t __bytes)
  {
    size_t __i = ((__bytes + (size_t)_S_align - 1) / (size_t)_S_align - 1);
    return _S_free_list + __i;
  }

  template <int Alignment>
    inline __mutex&
    pool_allocator_base<Alignment>::_M_get_mutex()
    { return palloc_init_mutex; }

// Allocate memory in large chunks in order to avoid fragmenting the
// heap too much.  Assume that __n is properly aligned.  We hold the
// allocation lock.
  template <int Alignment>
    inline char*
    pool_allocator_base<Alignment>::_M_allocate_chunk(size_t __n, int& __nobjs)
  {
    char* __result;
    size_t __total_bytes = __n * __nobjs;
    size_t __bytes_left = _S_end_free - _S_start_free;

    if (__bytes_left >= __total_bytes)
    {
      __result = _S_start_free;
      _S_start_free += __total_bytes;
      return __result ;
    }
    else if (__bytes_left >= __n)
    {
      __nobjs = (int)(__bytes_left / __n);
      __total_bytes = __n * __nobjs;
      __result = _S_start_free;
      _S_start_free += __total_bytes;
      return __result;
    }
    else
    {
// Try to make use of the left-over piece.
      if (__bytes_left > 0)
      {
        _Obj* volatile* __free_list = _M_get_free_list(__bytes_left);
        ((_Obj*)(void*)_S_start_free)->_M_free_list_link = *__free_list;
        *__free_list = (_Obj*)(void*)_S_start_free;
      }

      size_t __bytes_to_get = (2 * __total_bytes
        + _M_round_up(_S_heap_size >> 4));
      try
      {
        _S_start_free = static_cast<char*>(_mm_malloc(__bytes_to_get,Alignment));
      }
      catch (...)
      {
// Try to make do with what we have.  That can't hurt.  We
// do not try smaller requests, since that tends to result
// in disaster on multi-process machines.
        size_t __i = __n;
        for (; __i <= (size_t) _S_max_bytes; __i += (size_t) _S_align)
        {
          _Obj* volatile* __free_list = _M_get_free_list(__i);
          _Obj* __p = *__free_list;
          if (__p != 0)
          {
            *__free_list = __p->_M_free_list_link;
            _S_start_free = (char*)__p;
            _S_end_free = _S_start_free + __i;
            return _M_allocate_chunk(__n, __nobjs);
// Any leftover piece will eventually make it to the
// right free list.
          }
        }
// What we have wasn't enough.  Rethrow.
        _S_start_free = _S_end_free = 0;          // We have no chunk.
        __throw_exception_again;
      }
      _S_heap_size += __bytes_to_get;
      _S_end_free = _S_start_free + __bytes_to_get;
      return _M_allocate_chunk(__n, __nobjs);
    }
  }

// Returns an object of size __n, and optionally adds to "size
// __n"'s free list.  We assume that __n is properly aligned.  We
// hold the allocation lock.
  template <int Alignment>
    inline void*
    pool_allocator_base<Alignment>::_M_refill(size_t __n)
  {
    int __nobjs = 20;
    char* __chunk = _M_allocate_chunk(__n, __nobjs);
    _Obj* volatile* __free_list;
    _Obj* __result;
    _Obj* __current_obj;
    _Obj* __next_obj;

    if (__nobjs == 1)
      return __chunk;
    __free_list = _M_get_free_list(__n);

// Build free list in chunk.
    __result = (_Obj*)(void*)__chunk;
    *__free_list = __next_obj = (_Obj*)(void*)(__chunk + __n);
    for (int __i = 1; ; __i++)
    {
      __current_obj = __next_obj;
      __next_obj = (_Obj*)(void*)((char*)__next_obj + __n);
      if (__nobjs - 1 == __i)
      {
        __current_obj->_M_free_list_link = 0;
        break;
      }
      else
        __current_obj->_M_free_list_link = __next_obj;
    }
    return __result;
  }

  template <int Alignment>
    typename pool_allocator_base<Alignment>::_Obj* volatile pool_allocator_base<Alignment>::_S_free_list[_S_free_list_size];

  template <int Alignment>
    char* pool_allocator_base<Alignment>::_S_start_free = 0;

  template <int Alignment>
    char* pool_allocator_base<Alignment>::_S_end_free = 0;

  template <int Alignment>
    size_t pool_allocator_base<Alignment>::_S_heap_size = 0;

/// Class pool_allocator.
  template<typename _Tp, int Alignment = 8>
    class pool_allocator : private pool_allocator_base<Alignment>
  {
    private:
      static _Atomic_word     _S_force_new;

    public:
      typedef size_t                          size_type;
      typedef ptrdiff_t                       difference_type;
      typedef _Tp*                            pointer;
      typedef const _Tp*                      const_pointer;
      typedef _Tp&                            reference;
      typedef const _Tp&                      const_reference;
      typedef _Tp                             value_type;
      typedef pool_allocator_base<Alignment>  ancestor;
      typedef typename ancestor::_Obj         _Obj;

      template<typename _Tp1>
        struct rebind
        { typedef pool_allocator<_Tp1,Alignment> other; };

      pool_allocator() throw() { }

      pool_allocator(const pool_allocator&) throw() { }

      template<typename _Tp1>
        pool_allocator(const pool_allocator<_Tp1,Alignment>&) throw() { }

      ~pool_allocator() throw() { }

      pointer
        address(reference __x) const { return &__x; }

      const_pointer
        address(const_reference __x) const { return &__x; }

      size_type
        max_size() const throw()
        { return size_t(-1) / sizeof(_Tp); }

// _GLIBCXX_RESOLVE_LIB_DEFECTS
// 402. wrong new expression in [some_] allocator::construct
      void
        construct(pointer __p, const _Tp& __val)
        {
          ::new(__p) value_type(__val);
        }

      void
        destroy(pointer __p) { __p->~_Tp(); }

      pointer
        allocate(size_type __n, const void* = 0);

      void
        deallocate(pointer __p, size_type __n);
  };

  template<typename _Tp, int Alignment>
    inline bool
    operator==(const pool_allocator<_Tp,Alignment>&, const pool_allocator<_Tp,Alignment>&)
    { return true; }

  template<typename _Tp, int Alignment>
    inline bool
    operator!=(const pool_allocator<_Tp,Alignment>&, const pool_allocator<_Tp,Alignment>&)
    { return false; }

  template<typename _Tp, int Alignment>
    _Atomic_word
    pool_allocator<_Tp,Alignment>::_S_force_new;

  template<typename _Tp, int Alignment>
    inline _Tp*
    pool_allocator<_Tp,Alignment>::allocate(size_type __n, const void*)
  {
    pointer __ret = 0;
    if (__builtin_expect(__n != 0, true))
    {
      if (__builtin_expect(__n > this->max_size(), false))
        std::__throw_bad_alloc();

// If there is a race through here, assume answer from getenv
// will resolve in same direction.  Inspired by techniques
// to efficiently support threading found in basic_string.h.
      if (_S_force_new == 0)
      {
        if (std::getenv("GLIBCXX_FORCE_NEW"))
          __atomic_add_dispatch(&_S_force_new, 1);
        else
          __atomic_add_dispatch(&_S_force_new, -1);
      }

      const size_t __bytes = __n * sizeof(_Tp);
      if (__bytes > size_t(ancestor::_S_max_bytes) || _S_force_new == 1)
        __ret = static_cast<_Tp*>(_mm_malloc(__bytes,Alignment));
      else
      {
        _Obj* volatile* __free_list = ancestor::_M_get_free_list(__bytes);

        __scoped_lock sentry(ancestor::_M_get_mutex());
        _Obj* __restrict__ __result = *__free_list;
        if (__builtin_expect(__result == 0, 0))
          __ret = static_cast<_Tp*>(_M_refill(ancestor::_M_round_up(__bytes)));
        else
        {
          *__free_list = __result->_M_free_list_link;
          __ret = reinterpret_cast<_Tp*>(__result);
        }
        if (__builtin_expect(__ret == 0, 0))
          std::__throw_bad_alloc();
      }
    }
    return __ret;
  }

  template<typename _Tp, int Alignment>
    inline void
    pool_allocator<_Tp,Alignment>::deallocate(pointer __p, size_type __n)
  {
    if (__builtin_expect(__n != 0 && __p != 0, true))
    {
      const size_t __bytes = __n * sizeof(_Tp);
      if (__bytes > static_cast<size_t>(ancestor::_S_max_bytes) || _S_force_new == 1)
        _mm_free(__p);
      else
      {
        _Obj* volatile* __free_list = ancestor::_M_get_free_list(__bytes);
        _Obj* __q = reinterpret_cast<_Obj*>(__p);

        __scoped_lock sentry(ancestor::_M_get_mutex());
        __q ->_M_free_list_link = *__free_list;
        *__free_list = __q;
      }
    }
  }
}
#endif
