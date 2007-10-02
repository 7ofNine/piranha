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
      };
    public:
      simd_array(int n)
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            private_sub_[i].s[j]=n;
          }
        }
        for (i=0;i<rem;++i)
        {
          private_sub_[m128n-1].s[i]=n;
        }
        std::cout << "Features: " << (m128n) << ',' << (rem) << '\n';
      }
      simd_array(const simd_array &s)
      {
std::cout << "Copy ctor start\n";
        for (usint i=0;i<m128n;++i)
        {
          private_sub_[i].m=s.private_sub_[i].m;
        }
std::cout << "Copy ctor end\n";
      }
      simd_array &operator=(const simd_array &s)
      {
std::cout << "Assignment operator\n";
        if (this != &s)
        {
          for (usint i=0;i<m128n;++i)
          {
            private_sub_[i].m=s.private_sub_[i].m;
          }
        }
        return *this;
      }
      simd_array &operator+=(const simd_array &s)
      {
        usint i;
        for (i=0;i<m128n;++i)
        {
          private_sub_[i].m=_mm_add_epi16(private_sub_[i].m,s.private_sub_[i].m);
        }
        return *this;
      }
      void print(std::ostream &str) const
      {
        usint i,j;
        for (i=0;i<m128n-1;++i)
        {
          for (j=0;j<8;++j)
          {
            str << private_sub_[i].s[j] << '\t';
          }
        }
        for (i=0;i<rem;++i)
        {
          str << private_sub_[m128n-1].s[i] << '\t';
        }
      }
    private:
      static const sint                       m128n = Dim/8+1;
      static const sint                       dimenstion = Dim;
      static const sint                       rem = Dim%8;
      __attribute__ ((aligned (16))) subarray private_sub_[m128n];
  };
}

namespace std
{
  template <int Dim>
    ostream &operator<<(ostream &str, const simd_array<Dim> &s)
  {
    s.print(str);
    return str;
  }
}

#endif
