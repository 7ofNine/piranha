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

#include <boost/integer_traits.hpp>
#include <memory>
#include <string>
#include <vector>

#include "../base_classes/int_array.h"
#include "trig_array_commons.h"

namespace piranha
{
  template <int Bits, int Pos, class Allocator = std::allocator<char> >
  /// Trigonometric array, dynamically sized version.
  /**
   * It wraps a piranha::int_array with signed integer sized Bits, and adds the
   * capabilities needed for trigonometric manipulation.
   */
    class trig_array: public int_array<Bits,true,Allocator>,
    public trig_array_commons<trig_array<Bits,Pos>,Pos>
  {
      typedef trig_array_commons<trig_array,Pos> trig_commons;
      typedef int_array<Bits,true,Allocator> ancestor;
    public:
      typedef typename ancestor::value_type value_type;
      typedef typename ancestor::size_type size_type;
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
      template <class Series>
        bool checkup(const Series &s) const
      {
        switch (s.trig_width() != ancestor::size())
        {
          case true:
            std::cout << "Size mismatch in trig_array." << std::endl;
            return false;
          default:
            return true;
        }
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
      template <class ResultType>
        void trigmult(const trig_array &t2, ResultType &ret) const
      // NOTE: we are not using here a general version of vector addition/subtraction
      // because this way we can do two operations (+ and -) every cycle. This is a performance
      // critical part, so the optimization should be worth the hassle.
      {
        const size_type max_w=ancestor::size(), min_w=t2.size();
        // Assert widths, *this should always come from a regular Poisson series, and its width should hence be
        // already adjusted my merge_args in multiplication routines.
        p_assert(max_w >= min_w);
        p_assert(ret.template get<0>().size() == max_w);
        p_assert(ret.template get<1>().size() == max_w);
        size_type i;
        for (i=0;i<min_w;++i)
        {
          ret.template get<0>()[i]=(*this)[i]-t2[i];
          ret.template get<1>()[i]=(*this)[i]+t2[i];
        }
        for (;i<max_w;++i)
        {
          ret.template get<0>()[i]=(*this)[i];
          ret.template get<1>()[i]=(*this)[i];
        }
      }
      /// Does trig_array needs padding to be inserted in series with trig_width equal to n?
      template <class ArgsTuple>
        bool needs_padding(const ArgsTuple &args_tuple) const
      {
        return (ancestor::size() < args_tuple.template get<Pos>().size());
      }
      /// Does is trig_array insertable in series with trig_width equal to n?
      template <class ArgsTuple>
        bool is_insertable(const ArgsTuple &args_tuple) const
      {
        return (ancestor::size() <= args_tuple.template get<Pos>().size());
      }
      /// Pad right.
      template <class ArgsTuple>
        void pad_right(const ArgsTuple &args_tuple)
      {
        const size_t n=args_tuple.template get<Pos>().size();
        hard_assert(n >= ancestor::size());
        ancestor::resize(n);
      }
      // End INTERFACE definition.
      //-------------------------------------------------------
  };
}
#endif
