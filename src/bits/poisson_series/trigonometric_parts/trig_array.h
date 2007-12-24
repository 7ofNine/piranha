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

#ifndef PIRANHA_TRIG_ARRAY_H
#define PIRANHA_TRIG_ARRAY_H

#include <boost/integer_traits.hpp>

#include "../base_classes/base_trig_array.h"
#include "../../base_classes/int_array.h"

namespace piranha
{
  template <int Bits>
/// Trigonometric array, dynamically sized version.
    class trig_array: public base_trig_array<Bits,trig_array<Bits> >
  {
      typedef base_trig_array<Bits,trig_array> ancestor;
      typedef int_array<Bits,true> container_type;
      friend class base_trig_array<Bits,trig_array>;
    public:
      typedef typename ancestor::value_type value_type;
    public:
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_array():ancestor::base_trig_array() {}
/// Copy ctor.
      trig_array(const trig_array &t):ancestor::base_trig_array(t),private_container_(t.private_container_) {}
/// Ctor from piranha::deque_string.
      trig_array(const deque_string &sd):ancestor::base_trig_array()
      {
// TODO: check here that we are not loading too many multipliers, outside trig_size_t range.
// TODO: do it everywhere!
        const trig_size_t w=sd.size();
        if (w == 0)
        {
          std::cout << "Warning: constructing empty trig_array." << std::endl;
          std::abort();
          return;
        }
// Now we know  w >= 1.
        private_container_.resize(w-1);
        for (trig_size_t i=0;i < w-1;++i)
        {
          private_container_[i]=utils::lexical_converter<value_type>(sd[i]);
        }
// Take care of flavour.
        if (*sd.back().c_str() == 's')
        {
          ancestor::flavour()=false;
        }
      }
      ~trig_array() {}
      void pad_right(const size_t &n)
      {
        p_assert(n >= private_container_.size());
        private_container_.resize(n);
      }
// Probing.
/// Data footprint.
/**
 * Returns the memory occupied by the data members.
 */
      size_t data_footprint() const
      {
        return (g_width()*sizeof(value_type));
      }
      template <class Series>
        bool checkup(const Series &s) const
      {
        if (s.trig_width() != g_width())
        {
          std::cout << "Size mismatch in trig_array." << std::endl;
          return false;
        }
        return true;
      }
      bool needs_padding(const size_t &n) const
      {
        return (g_width() < n);
      }
      bool is_insertable(const size_t &n) const
      {
        return (g_width() <= n);
      }
// FIXME: introduce size_type from int_array here.
      static const size_t max_size = boost::integer_traits<size_t>::const_max;
      bool operator==(const trig_array &t2) const
      {
        return (ancestor::flavour() == t2.flavour() and private_container_ == t2.private_container_);
      }
      bool operator<(const trig_array &t2) const
      {
        return ancestor::less_than(t2);
      }
      size_t hasher() const
      {
        size_t retval(private_container_.hasher());
        boost::hash_combine(retval,ancestor::flavour());
        return retval;
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
      void trigmult(const trig_array &t2, trig_array &ret1, trig_array &ret2) const
// NOTE: we are not using here a general version of vector addition/subtraction
// because this way we can do two operations (+ and -) every cycle. This is a performance
// critical part, so the optimization should be worth the hassle.
      {
        const trig_size_t max_w=g_width(), min_w=t2.g_width();
// Assert widths, *this should always come from a regular Poisson series, and its width should hence be
// already adjusted my merge_args in multiplication routines.
        p_assert(max_w >= min_w);
        p_assert(ret1.g_width() == max_w);
        p_assert(ret2.g_width() == max_w);
        trig_size_t i;
        for (i=0;i<min_w;++i)
        {
          ret1.private_container_[i]=private_container_[i]-t2.private_container_[i];
          ret2.private_container_[i]=private_container_[i]+t2.private_container_[i];
        }
        for (;i<max_w;++i)
        {
          ret1.private_container_[i]=private_container_[i];
          ret2.private_container_[i]=private_container_[i];
        }
      }
      trig_array &operator*=(const int &n)
      {
        ancestor::mult_by_int(n);
        return *this;
      }
// End INTERFACE definition.
//-------------------------------------------------------
    private:
      trig_size_t g_width() const
      {
        return private_container_.size();
      }
      const value_type *g_container() const
      {
        return &(private_container_[0]);
      }
      value_type *s_container()
      {
        return &(private_container_[0]);
      }
// Data members.
    private:
      container_type        private_container_;
  };
}
#endif
