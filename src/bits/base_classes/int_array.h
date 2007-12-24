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

#ifndef PIRANHA_INT_ARRAY_H
#define PIRANHA_INT_ARRAY_H

#include <boost/integer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <memory> // For std::allocator.

#include "../p_assert.h"

namespace piranha
{
// Used below to define the integer type.
  template <int Bits, bool Signed>
    struct integer_type_chooser
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    typedef typename boost::int_t<Bits>::fast type;
    typedef max_fast_int max_fast_type;
  };

// Specialization for unsigned integers.
  template <int Bits>
    struct integer_type_chooser<Bits,false>
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    typedef typename boost::uint_t<Bits>::fast type;
    typedef max_fast_uint max_fast_type;
  };

/// Meta-programmed functor for the calculation of base-2 logarithm.
/**
 * Result is retrieved through the lg::value const member function.
 */
  template <int N, int Cur = 0>
    struct lg
  {
    BOOST_STATIC_ASSERT(N > 0 and (N % 2) == 0);
    static const size_t value = lg<(N >> 1),Cur+1>::value;
  };

  template <int Cur>
    struct lg<1,Cur>
  {
    static const size_t value = Cur;
  };

/// Dynamically-sized integer array.
/**
 * Parametrized to an integer sized Bits, which can be Signed or not.
 */
  template <int Bits, bool Signed, class Allocator = std::allocator<char> >
    class int_array
  {
      typedef typename integer_type_chooser<Bits,Signed>::type value_type;
      typedef typename integer_type_chooser<Bits,Signed>::max_fast_type max_fast_type;
      typedef uint16 size_type;
      typedef typename Allocator::template rebind<value_type>::other allocator_type;
      BOOST_STATIC_ASSERT(sizeof(max_fast_type) % sizeof(value_type) == 0);
    public:
/// Cast argument to int_array::max_fast_type pointer.
#define max_cast(arg) ((max_fast_type *)(arg))
/// Default ctor.
/**
 * Constructs empty array.
 */
      int_array():m_size(0),m_pack_size(0),m_ptr(allocator.allocate(0)) {}
/// Copy ctor.
      int_array(const int_array &v):m_size(v.m_size),m_pack_size(v.m_pack_size),m_ptr(allocator.allocate(m_size))
      {
        packed_copy(m_ptr,v.m_ptr,m_size,m_pack_size);
      }
/// Ctor from size.
/**
 * Initialises to 0 all the elements.
 */
      int_array(const size_type &s):m_size(s),m_pack_size(s >> pack_shift),m_ptr(allocator.allocate(m_size))
      {
        size_type i;
        for (i=0;i < m_pack_size;++i)
        {
          max_cast(m_ptr)[i]=0;
        }
        for(i = i << pack_shift;i < m_size;++i)
        {
          m_ptr[i]=0;
        }
      }
/// Dtor.
      ~int_array() {allocator.deallocate(m_ptr,m_size);}
      const value_type &operator[](const size_t &n) const {return m_ptr[n];}
      value_type &operator[](const size_t &n) {return m_ptr[n];}
      size_t size() const {return m_size;}
      void resize(const size_type &new_size)
      {
        value_type *new_ptr = alloc_if_size_differs(new_size);
        if (new_ptr == m_ptr)
        {
          return;
        }
        const size_type new_pack_size = (new_size >> pack_shift);
// Copy to the minimum of the new sizes.
        packed_copy(new_ptr,m_ptr,std::min(m_size,new_size),std::min(m_pack_size,new_pack_size));
// Zero the remaining elements, if any.
        for (size_type i=m_size;i < new_size;++i)
        {
          new_ptr[i]=0;
        }
// Destroy old pointer and assign new data members.
        allocator.deallocate(m_ptr,m_size);
        m_ptr = new_ptr;
        m_size = new_size;
        m_pack_size = new_pack_size;
      }
/// Assignment operator.
      int_array &operator=(const int_array &v)
      {
        value_type *new_ptr = alloc_if_size_differs(v.m_size);
        if (new_ptr != m_ptr)
        {
          allocator.deallocate(m_ptr,m_size);
          m_ptr = new_ptr;
          m_size = v.m_size;
          m_pack_size = v.m_pack_size;
        }
        packed_copy(m_ptr,v.m_ptr,m_size,m_pack_size);
        return *this;
      }
/// Hash value.
      size_t hasher() const
      {
        size_t retval=0;
        size_type i;
        for (i=0;i < m_pack_size;++i)
        {
          boost::hash_combine(retval,max_cast(m_ptr)[i]);
        }
        for(i = i << pack_shift;i < m_size;++i)
        {
          boost::hash_combine(retval,m_ptr[i]);
        }
        return retval;
      }
/// Equality test.
      bool operator==(const int_array &v) const
      {
        switch (m_size == v.size())
        {
          case true:
            {
              size_type i;
              for (i=0;i < m_pack_size;++i)
              {
                if (max_cast(m_ptr)[i] != max_cast(v.m_ptr)[i])
                {
                  return false;
                }
              }
              for(i = i << pack_shift;i < m_size;++i)
              {
                if (m_ptr[i] != v[i])
                {
                  return false;
                }
              }
            }
            return true;
// TODO: experiment here with case false?
          default:
            return false;
        }
      }
/// Test for zero elements.
      bool is_zero() const
      {
        size_t i;
        for (i=0;i < m_pack_size;++i)
        {
          if (max_cast(m_ptr)[i] != 0)
          {
            return false;
          }
        }
        for(i = i << pack_shift;i < m_size;++i)
        {
          if (m_ptr[i] != 0)
          {
            return false;
          }
        }
        return true;
      }
    private:
      void packed_copy(value_type *new_ptr, const value_type *old_ptr, const size_type &size,
        const size_type &pack_size)
      {
        size_type i;
        for (i=0;i < pack_size;++i)
        {
          max_cast(new_ptr)[i]=max_cast(old_ptr)[i];
        }
        for(i = i << pack_shift;i < size;++i)
        {
          new_ptr[i]=old_ptr[i];
        }
      }
// Allocate if new_size is different from current size.
      value_type *alloc_if_size_differs(const size_type &new_size)
      {
        switch (m_size == new_size)
        {
          case true:
            return m_ptr;
          default:
            return allocator.allocate(new_size);
        }
      }
    private:
// Data members.
/// Size of the array.
/**
 * Equal to the number of elements contained by the array.
 */
      size_type               m_size;
/// Packed size of the array.
/**
 * Defined by the integer division int_array::m_size / int_array::pack_mult.
 */
      size_type               m_pack_size;
/// Pointer to the first value of the array.
      value_type              *m_ptr;
/// Array allocator.
      static allocator_type   allocator;
/// Pack multiplier.
/**
 * Defined by the integer division between the number of bits of int_array::max_fast_type and the number
 * of bits of int_array::value_type.
 */
      static const size_type  pack_mult = sizeof(max_fast_type)/sizeof(value_type);
/// Pack shifting.
/**
 * Defined as the base-2 logarithm of int_array::pack_mult. If int_array::pack_mult is not a power of two,
 * a compilation error is produced.
 */
      static const size_type  pack_shift = lg<pack_mult>::value;
#undef max_cast
  };

  template <int Bits, bool Signed, class Allocator>
    typename int_array<Bits,Signed,Allocator>::allocator_type int_array<Bits,Signed,Allocator>::allocator;
};

#endif
