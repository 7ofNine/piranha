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

#ifndef PIRANHA_TRIG_FIXED_ARRAY_H
#define PIRANHA_TRIG_FIXED_ARRAY_H

#include <boost/algorithm/minmax.hpp>
#include <cstring>
#include <iostream>

#include "../base_classes/trig_array_commons.h"
#include "../../common_typedefs.h"
#include "../../packed_int_array.h"
#include "../../type_traits/is_resizable.h"
#include "../../type_traits/false_type.h"
#include "../../utils.h"

namespace piranha
{
  template <int N>
    class tfa_unrollers
  {
    public:
      template <class T>
        static void mult(const T *end_array1, const T *end_array2, T *ret_array1, T *ret_array2)
      {
        ret_array1[-N]=end_array1[-N]-end_array2[-N];
        ret_array2[-N]=end_array1[-N]+end_array2[-N];
        tfa_unrollers<N-1>::mult(end_array1,end_array2,ret_array1,ret_array2);
      }
      template <class T>
        static void invert(T *end_array)
      {
        end_array[-N]=-end_array[-N];
        tfa_unrollers<N-1>::invert(end_array);
      }
  };

  template <>
    class tfa_unrollers<1>
  {
    public:
      template <class T>
        static void mult(const T *end_array1, const T *end_array2, T *ret_array1, T *ret_array2)
      {
        ret_array1[-1]=end_array1[-1]-end_array2[-1];
        ret_array2[-1]=end_array1[-1]+end_array2[-1];
      }
      template <class T>
        static void invert(T *end_array)
      {
        end_array[-1]=-end_array[-1];
      }
  };

/// Trigonometric array, fixed size version.
  template <int Dim, int Bits>
    class trig_fixed_array: public packed_int_array<Dim,Bits>,
      public trig_array_commons<trig_fixed_array<Dim,Bits> >
  {
      typedef packed_int_array<Dim,Bits> ancestor;
      typedef trig_array_commons<trig_fixed_array<Dim,Bits> > trig_commons;
    public:
      typedef typename ancestor::value_type value_type;
      typedef typename ancestor::size_type size_type;
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_fixed_array():ancestor::packed_int_array() {}
/// Ctor from piranha::deque_string.
      trig_fixed_array(const deque_string &sd):ancestor::packed_int_array(),trig_commons::trig_array_commons(sd) {}
      static const size_t &size() {return ancestor::max_size;}
// TODO: check if it is better this or the base version.
      /*void invert_sign()
      {
        tfa_unrollers<dimension>::invert(&private_container_[0]+dimension);
      }*/
// Probing.
/// Data footprint.
/**
 * Returns the memory occupied by the data members.
 */
      size_t data_footprint() const
      {
        return sizeof(trig_fixed_array);
      }
      template <class Series>
        bool checkup(const Series &s) const
      {
        return s.cf_width() <= Dim;
      }
      bool needs_padding(const size_t &n) const
      {
        p_assert(n <= Dim);
// Disable compiler warning when asserts are disabled.
        (void)n;
        return false;
      }
      bool is_insertable(const size_t &n) const
      {
        return (Dim >= n);
      }
// All multipliers are zero.
//TODO: use packing here and place this method in packed int array.
      bool is_zero() const
      {
        for (size_type i=0;i<Dim;++i)
        {
          if ((*this)[i]!=0)
          {
            return false;
          }
        }
        return true;
      }
// Math.
/// Multiplication.
/**
 * Multiplication of two trigonometric functions using Werner's formulas, i.e.
 * \f[
 * C\cos\alpha\cdot\cos\beta=
 * \frac{C}{2} \cos \left( \alpha - \beta \right) + \frac{C}{2} \cos \left( \alpha + \beta \right)
 * \f]
 * and the likes. Notice that in the first return value always goes the \f$ \alpha - \beta \f$ term
 * and in the second one always goes \f$ \alpha + \beta \f$ one.
 * Please also note that no assumptions are made with respect to return values' content (e.g., it is not guaranteed
 * that return values are empty).
 * @param[in] t2 factor.
 * @param[out] ret1 first return value.
 * @param[out] ret2 second return value.
 */
      void trigmult(const trig_fixed_array &t2, trig_fixed_array &ret1, trig_fixed_array &ret2) const
      {
#ifdef _PIRANHA_SSE2
        ancestor::trigmult(*this,t2,ret1,ret2);
#else
        tfa_unrollers<Dim>::mult(&(*this)[0]+Dim,&t2[0]+Dim,&ret1[0]+Dim,&ret2[0]+Dim);
#endif
      }
// End INTERFACE definition.
//-------------------------------------------------------
  };

/// Resizable type-traits specialization for piranha::trig_fixed_array.
  template <>
    template <int Dim, int Bits>
    struct is_resizable<trig_fixed_array<Dim,Bits> >:public false_type {};
}

#endif
