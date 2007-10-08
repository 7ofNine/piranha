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

#ifndef PIRANHA_PACKED_INT_ARRAY32_H
#define PIRANHA_PACKED_INT_ARRAY32_H

#include <boost/integer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>

#include "common_typedefs.h"

namespace piranha
{
// TODO: check the real impact of these unrollers.
  template <int N>
    inline void hash_unroller(size_t &seed, const int16 *end_array)
  {
    boost::hash_combine(seed,end_array[-N]);
    hash_unroller<N-1>(seed,end_array);
  }

  template <>
    inline void hash_unroller<1>(size_t &seed, const int16 *end_array)
  {
    boost::hash_combine(seed,end_array[-1]);
  }

  template <int Dim, int Bits>
    class packed_int_array
  {
      BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
      BOOST_STATIC_ASSERT(Dim > 0 and Dim < 100);
    public:
      typedef typename boost::int_t<Bits>::fast value_type;
      packed_int_array()
        {}
      void hasher(size_t &seed) const
      {
        uint8 i;
        for (i=0;i<size32;++i)
        {
//std::cout << "Entered 32\n";
          boost::hash_combine(seed,*((const int32 *)(const void *)(private_container_)+i));
        }
        for (i=(size32<<1);i<size16;++i)
        {
//std::cout << "Entered 16\n";
          boost::hash_combine(seed,*((const int16 *)(const void *)(private_container_)+i));
        }
        for (i=(size16<<1);i<size8;++i)
        {
//std::cout << "Entered 8\n";
          boost::hash_combine(seed,*((const int8 *)(const void *)(private_container_)+i));
        }
      }
      bool operator==(const packed_int_array &p) const
      {
        uint8 i;
        for (i=0;i<size32;++i)
        {
          if (*((const int32 *)(const void *)(private_container_)+i) !=
            *((const int32 *)(const void *)(p.private_container_)+i))
          {
            return false;
          }
        }
        for (i=(size32<<1);i<size16;++i)
        {
          if (*((const int16 *)(const void *)(private_container_)+i) !=
            *((const int16 *)(const void *)(p.private_container_)+i))
          {
            return false;
          }
        }
        for (i=(size16<<1);i<size8;++i)
        {
          if (*((const int8 *)(const void *)(private_container_)+i) !=
            *((const int8 *)(const void *)(p.private_container_)+i))
          {
            return false;
          }
        }
        return true;
      }
      const value_type &operator[](uint8 n) const
      {
        return private_container_[n];
      }
      value_type &operator[](uint8 n)
      {
        return private_container_[n];
      }
    private:
      value_type          private_container_[Dim];
      static const uint8  size32 = (uint8)((Dim*Bits)/32);
      static const uint8  size16 = (uint8)((Dim*Bits)/16);
      static const uint8  size8 = (uint8)((Dim*Bits)/8);
  };

  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size32;

  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size16;

  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size8;
}

#endif
