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

#ifndef PIRANHA_SIMD_ARRAY_H
#define PIRANHA_SIMD_ARRAY_H

#include <boost/static_assert.hpp>
#include <iostream>
#include <emmintrin.h>

#include "common_typedefs.h"

// FIXME: here we must pay attention to the fact that mult_t != sint, and on other arch this could
// be an issue. Maybe we should use the types provided by the compiler? uint32, uint16, etc.

namespace piranha
{
  template <int Dim>
    class simd_array
  {
      BOOST_STATIC_ASSERT(Dim > 0);
      BOOST_STATIC_ASSERT(Dim < 100);
      union subarray
      {
        __m128i m;
        sint    s[8];
        int     i[4];
      };
    public:
/// Default constructor.
/**
 * Values will be uninitialized.
 **/
      simd_array()
        {}
/// Copy constructor.
      simd_array(const simd_array &s)
      {
//std::cout << "Copy ctor start\n";
        assignment(s);
//std::cout << "Copy ctor end\n";
      }
/// Destructor.
      ~simd_array()
        {}
// These operators are used from piranha::base_trig_array in non-speed-critical operations (I/O, for instance).
      const sint &operator[](const usint &n) const
      {
        return private_sub_[n/8].s[n%8];
      }
      sint &operator[](const usint &n)
      {
        return private_sub_[n/8].s[n%8];
      }
      simd_array &operator=(const simd_array &s)
      {
//std::cout << "Assignment operator\n";
        if (this != &s)
        {
          assignment(s);
        }
        return *this;
      }
      static void add(const simd_array &s1, const simd_array &s2, simd_array &result)
      {
        usint i;
        for (i=0;i<m128n;++i)
        {
          result.private_sub_[i].m=_mm_add_epi16(s1.private_sub_[i].m,s2.private_sub_[i].m);
        }
      }
      static void sub(const simd_array &s1, const simd_array &s2, simd_array &result)
      {
        usint i;
        for (i=0;i<m128n;++i)
        {
          result.private_sub_[i].m=_mm_sub_epi16(s1.private_sub_[i].m,s2.private_sub_[i].m);
        }
      }
      bool equality_test(const simd_array &s) const
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<4;++j)
          {
            if (private_sub_[i].i[j] != s.private_sub_[i].i[j])
            {
              return false;
            }
          }
        }
        for (i=0;i<rem;++i)
        {
          if (private_sub_[m128n-1].s[i] != s.private_sub_[m128n-1].s[i])
          {
            return false;
          }
        }
        return true;
      }
      bool less_than(const simd_array &s) const
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            if (private_sub_[i].s[j] < s.private_sub_[i].s[j])
            {
              return true;
            }
            else if (private_sub_[i].s[j] > s.private_sub_[i].s[j])
            {
              return false;
            }
          }
        }
        for (i=0;i<rem;++i)
        {
          if (private_sub_[m128n-1].s[i] < s.private_sub_[m128n-1].s[i])
          {
            return true;
          }
          else if (private_sub_[m128n-1].s[i] < s.private_sub_[m128n-1].s[i])
          {
            return false;
          }
        }
        return false;
      }
      void hasher(size_t &seed) const
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<4;++j)
          {
            boost::hash_combine(seed,private_sub_[i].i[j]);
          }
        }
        for (i=0;i<rem;++i)
        {
          boost::hash_combine(seed,private_sub_[m128n-1].s[i]);
        }
      }
      simd_array &operator*=(const sint &n)
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            private_sub_[i].s[j]*=n;
          }
        }
        for (i=0;i<rem;++i)
        {
          private_sub_[m128n-1].s[i]*=n;
        }
        return *this;
      }
      short int sign() const
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            if (private_sub_[i].s[j] > 0)
            {
              return 1;
            }
            if (private_sub_[i].s[j] < 0)
            {
              return -1;
            }
          }
        }
        for (i=0;i<rem;++i)
        {
          if (private_sub_[m128n-1].s[i] > 0)
          {
            return 1;
          }
          if (private_sub_[m128n-1].s[i] < 0)
          {
            return -1;
          }
        }
        return 1;
      }
      void invert_sign()
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            private_sub_[i].s[j]=-private_sub_[i].s[j];
          }
        }
        for (i=0;i<rem;++i)
        {
          private_sub_[m128n-1].s[i]=-private_sub_[m128n-1].s[i];
        }
      }
      bool is_zero() const
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            if (private_sub_[i].s[j] != 0)
            {
              return false;
            }
          }
        }
        for (i=0;i<rem;++i)
        {
          if (private_sub_[m128n-1].s[i] != 0)
          {
            return false;
          }
        }
        return true;
      }
    private:
      void assignment(const simd_array &s)
      {
        for (usint i=0;i<m128n;++i)
        {
          private_sub_[i].m=s.private_sub_[i].m;
        }
      }
    private:
      static const sint                       m128n = Dim/8+1;
      static const sint                       dimenstion = Dim;
      static const sint                       rem = Dim%8;
      __attribute__ ((aligned (16))) subarray private_sub_[m128n];
  };
}

#endif
