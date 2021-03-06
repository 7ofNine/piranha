/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
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

#ifndef PIRANHA_PACKED_INT_ARRAY_H
#define PIRANHA_PACKED_INT_ARRAY_H

#include <boost/integer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>

#include "common_typedefs.h"
#include "piranha_vectorization.h"

namespace piranha
{
  // TODO: check the real impact of these unrollers.
  //   template <int N>
  //     inline void hash_unroller(size_t &seed, const int16 *end_array)
  //   {
  //     boost::hash_combine(seed,end_array[-N]);
  //     hash_unroller<N-1>(seed,end_array);
  //   }
  //
  //   template <>
  //     inline void hash_unroller<1>(size_t &seed, const int16 *end_array)
  //   {
  //     boost::hash_combine(seed,end_array[-1]);
  //   }

  template <int Dim, int Bits> class packed_int_array;
  template <int Bits> class pia_helper;

#ifdef _PIRANHA_SSE2
  template <int Bits>
    class pia_helper
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    template <int Dim, int Bits_> friend class packed_int_array;
    // Default is 16bit operations.
    template <class T>
      static void trigmult(const T &t1, const T &t2, T &ret1, T &ret2)
    {
      typedef typename T::size_type size_type;
      for (size_type i=0;i< T::m128n;++i)
      {
        ret1.private_container_.m[i]=_mm_sub_epi16(t1.private_container_.m[i],t2.private_container_.m[i]);
        ret2.private_container_.m[i]=_mm_add_epi16(t1.private_container_.m[i],t2.private_container_.m[i]);
      }
    }
  };

  template <>
    class pia_helper<8>
  {
    template <int Dim, int Bits> friend class packed_int_array;
    // Default is 16bit operations.
    template <class T>
      static void trigmult(const T &t1, const T &t2, T &ret1, T &ret2)
    {
      typedef typename T::size_type size_type;
      for (size_type i=0;i < T::m128n;++i)
      {
        ret1.private_container_.m[i]=_mm_sub_epi8(t1.private_container_.m[i],t2.private_container_.m[i]);
        ret2.private_container_.m[i]=_mm_add_epi8(t1.private_container_.m[i],t2.private_container_.m[i]);
      }
    }
  };
#endif

  template <int Dim, int Bits>
    class packed_int_array
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    BOOST_STATIC_ASSERT(Dim > 0 and Dim < 100);
    public:
      typedef typename boost::int_t<Bits>::fast value_type;
      typedef uint8 size_type;
    private:
#ifdef _PIRANHA_SSE2
      template <int Bits_> friend class pia_helper;
      static const size_type m128n = (size_type)((Dim*Bits)/128+1);
#endif
      union container_type
      {
        value_type  v[Dim];
#ifdef _PIRANHA_SSE2
        // We want to have always at least 1 __m128i object to take advantage of SSE2.
        __attribute__ ((aligned (16))) __m128i     m[m128n];
#endif
      };
    public:
      packed_int_array():m_flavour(true)
      {
        for (size_type i=0;i<Dim;++i)
        {
          private_container_.v[i]=0;
        }
      }
      size_t hasher() const
      {
        size_t retval=0;
        size_type i;
#ifdef _PIRANHA_64BIT
        for (i=0;i<size64;++i)
        {
          boost::hash_combine(retval,*((const int64 *)(const void *)(private_container_.v)+i));
        }
        for (i=(size64<<1);i<size32;++i)
#else
          for (i=0;i<size32;++i)
#endif
      {
        boost::hash_combine(retval,*((const int32 *)(const void *)(private_container_.v)+i));
      }
      for (i=(size32<<1);i<size16;++i)
      {
        boost::hash_combine(retval,*((const int16 *)(const void *)(private_container_.v)+i));
      }
      for (i=(size16<<1);i<size8;++i)
      {
        boost::hash_combine(retval,*((const int8 *)(const void *)(private_container_.v)+i));
      }
      return retval;
    }
    bool equal_to(const packed_int_array &p) const
    {
      size_type i;
#ifdef _PIRANHA_64BIT
      for (i=0;i<size64;++i)
      {
        if (*((const int64 *)(const void *)(private_container_.v)+i) !=
          *((const int64 *)(const void *)(p.private_container_.v)+i))
        {
          return false;
        }
      }
      for (i=(size64<<1);i<size32;++i)
#else
        for (i=0;i<size32;++i)
#endif
      {
        if (*((const int32 *)(const void *)(private_container_.v)+i) !=
          *((const int32 *)(const void *)(p.private_container_.v)+i))
        {
          return false;
        }
      }
      for (i=(size32<<1);i<size16;++i)
      {
        if (*((const int16 *)(const void *)(private_container_.v)+i) !=
          *((const int16 *)(const void *)(p.private_container_.v)+i))
        {
          return false;
        }
      }
      for (i=(size16<<1);i<size8;++i)
      {
        if (*((const int8 *)(const void *)(private_container_.v)+i) !=
          *((const int8 *)(const void *)(p.private_container_.v)+i))
        {
          return false;
        }
      }
      return true;
    }
    const value_type &operator[](size_type n) const {return private_container_.v[n];}
    value_type &operator[](size_type n) {return private_container_.v[n];}
    const bool &flavour() const {return m_flavour;}
    bool &flavour() {return m_flavour;}
    static const size_t max_size = Dim;
    /// Resize.
    /**
     * Resizing is allowed as long as a size greater than Dim is requested.
     */
    void resize(const size_t &n)
    {
      hard_assert(n <= Dim);
      (void)n;
    }
#ifdef _PIRANHA_SSE2
    static void trigmult(const packed_int_array &t1, const packed_int_array &t2,
      packed_int_array &ret1, packed_int_array &ret2)
    {
      pia_helper<Bits>::trigmult(t1,t2,ret1,ret2);
    }
#endif
    private:
      bool                    m_flavour;
      container_type          private_container_;
      static const size_type  size64 = (uint8)((Dim*Bits)/64);
      static const size_type  size32 = (uint8)((Dim*Bits)/32);
      static const size_type  size16 = (uint8)((Dim*Bits)/16);
      static const size_type  size8 = (uint8)((Dim*Bits)/8);
  };

  // Static inits.
  template <int Dim, int Bits>
    const size_t packed_int_array<Dim,Bits>::max_size;

  template <int Dim, int Bits>
    const typename packed_int_array<Dim,Bits>::size_type packed_int_array<Dim,Bits>::size64;

  template <int Dim, int Bits>
    const typename packed_int_array<Dim,Bits>::size_type packed_int_array<Dim,Bits>::size32;

  template <int Dim, int Bits>
    const typename packed_int_array<Dim,Bits>::size_type packed_int_array<Dim,Bits>::size16;

  template <int Dim, int Bits>
    const typename packed_int_array<Dim,Bits>::size_type packed_int_array<Dim,Bits>::size8;

#ifdef _PIRANHA_SSE2
  template <int Dim, int Bits>
    const typename packed_int_array<Dim,Bits>::size_type packed_int_array<Dim,Bits>::m128n;
#endif
}
#endif
