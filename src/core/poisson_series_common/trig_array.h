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

#ifndef PIRANHA_TRIG_ARRAY_H
#define PIRANHA_TRIG_ARRAY_H

#include <memory> // For standard allocator.
#include <string>

#include "../base_classes/int_array.h"
#include "trig_array_commons.h"

#define __PIRANHA_TRIG_ARRAY_TP_DECL int Bits, int Pos, class Allocator
#define __PIRANHA_TRIG_ARRAY_TP Bits,Pos,Allocator

namespace piranha
{
  /// Trigonometric array, dynamically sized version.
  /**
   * It wraps a piranha::int_array with signed integer sized Bits, and adds the
   * capabilities needed for trigonometric manipulation.
   */
  template <__PIRANHA_TRIG_ARRAY_TP_DECL = std::allocator<char> >
    class trig_array:
    public int_array<Bits,Pos,true,Allocator,trig_array<__PIRANHA_TRIG_ARRAY_TP> >,
    public trig_array_commons<trig_array<__PIRANHA_TRIG_ARRAY_TP> >
  {
      friend class trig_array_commons<trig_array<__PIRANHA_TRIG_ARRAY_TP> >;
      typedef trig_array_commons<trig_array<__PIRANHA_TRIG_ARRAY_TP> > trig_commons;
      typedef int_array<Bits,Pos,true,Allocator,trig_array<__PIRANHA_TRIG_ARRAY_TP> > ancestor;
    public:
      typedef typename ancestor::value_type value_type;
      typedef typename ancestor::size_type size_type;
      typedef typename trig_commons::coded_type coded_type;
      // Start INTERFACE definition.
      //-------------------------------------------------------
      // Ctors.
      /// Default ctor.
      trig_array():ancestor::int_array() {}
      /// Ctor from string.
      template <class ArgsTuple>
        explicit trig_array(const std::string &s, const ArgsTuple &):ancestor::int_array(),trig_commons::trig_array_commons(s) {}
      // Probing.
      /// Data footprint.
      /**
       * Returns the memory occupied by the data members.
       */
      size_t data_footprint() const {return (ancestor::size()*sizeof(value_type));}
      // Math.
      /// Multiplication.
      /**
       * Used in poisson_series_term multiplication.
       * TODO: update docs below.
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
      void multiply(const trig_array &t2, trig_array &ret1, trig_array &ret2) const
      // NOTE: we are not using here a general version of vector addition/subtraction
      // because this way we can do two operations (+ and -) every cycle. This is a performance
      // critical part, so the optimization should be worth the hassle.
      {
        const size_type max_w=ancestor::size(), min_w=t2.size();
        // Assert widths, *this should always come from a regular Poisson series, and its width should hence be
        // already adjusted my merge_args in multiplication routines.
        p_assert(max_w >= min_w);
        // Adjust the width of retvals, if needed.
        ret1.resize(max_w);
        ret2.resize(max_w);
        p_assert(ret1.size() == max_w);
        p_assert(ret2.size() == max_w);
        size_type i;
        for (i=0; i < min_w; ++i)
        {
          ret1[i]=(*this)[i]-t2[i];
          ret2[i]=(*this)[i]+t2[i];
        }
        for (; i < max_w; ++i)
        {
          ret1[i]=(*this)[i];
          ret2[i]=(*this)[i];
        }
      }
      // End INTERFACE definition.
      //-------------------------------------------------------
  };

  // Representation for trigonometric array during series multiplication is a const pointer.
  template <__PIRANHA_TRIG_ARRAY_TP_DECL>
    struct series_mult_rep<trig_array<__PIRANHA_TRIG_ARRAY_TP> >
  {
    typedef trig_array<__PIRANHA_TRIG_ARRAY_TP> const * type;
    static const trig_array<__PIRANHA_TRIG_ARRAY_TP> &get(type p) {return *p;}
    static void assign(type &res, const trig_array<__PIRANHA_TRIG_ARRAY_TP> &source) {res=&source;}
  };
}

#undef __PIRANHA_TRIG_ARRAY_TP_DECL
#undef __PIRANHA_TRIG_ARRAY_TP

#endif
