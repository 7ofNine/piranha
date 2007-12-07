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

#include "../base_classes/base_trig_array.h"
#include "../bits/type_traits/is_resizable.h"
#include "../bits/type_traits/true_type.h"

namespace piranha
{
  template <int Bits>
/// Trigonometric array, dynamically sized version.
    class trig_array: public base_trig_array<Bits,trig_array<Bits> >
  {
      typedef base_trig_array<Bits,trig_array> ancestor;
    public:
      typedef typename ancestor::value_type value_type;
    private:
      typedef std::valarray<value_type> container_type;
      template <int Bits_, class Derived>
        friend class base_trig_array;
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
        if (w==0)
        {
          std::cout << "Warning: constructing empty trig_array." << std::endl;
          std::abort();
          return;
        }
// Now we know  w >= 1.
        private_container_.resize(w-1);
        for (trig_size_t i=0;i<w-1;++i)
        {
          private_container_[i]=utils::lexical_converter<value_type>(sd[i]);
        }
// Take care of flavour.
        if (*sd.back().c_str()=='s')
        {
          ancestor::s_flavour()=false;
        }
      }
      ~trig_array() {}
// Manip.
      void prepend_args(const size_t &n)
      {
        if (n>0)
        {
          container_type old_private_container_(private_container_);
          const trig_size_t old_w=old_private_container_.size();
          private_container_.resize(old_w+n);
          for (trig_size_t i=0;i<old_w;++i)
          {
            private_container_[i+n]=old_private_container_[i];
          }
        }
      }
      void append_args(const size_t &n)
      {
        increase_size(g_width()+n);
      }
      void increase_size(const size_t &n)
      {
        p_assert(n >= g_width());
        if (n > g_width())
        {
          container_type old_private_container_(private_container_);
          private_container_.resize(n);
          const trig_size_t old_w=old_private_container_.size();
          for (trig_size_t i=0;i<old_w;++i)
          {
            private_container_[i]=old_private_container_[i];
          }
        }
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
      bool checkup(const size_t &w) const
      {
        if (w != g_width())
        {
          std::cout << "Size mismatch in trig_array." << std::endl;
          return false;
        }
        return true;
      }
      bool operator==(const trig_array &t2) const
      {
        return ancestor::equality_test(t2);
      }
      bool operator<(const trig_array &t2) const
      {
        return ancestor::less_than(t2);
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
      trig_array &operator=(const trig_array &t2)
      {
        ancestor::assignment_operator(t2);
        return *this;
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
/// Assignment.
/**
 * After the assignment the data members of the two classes must be equal, both from a mathematical
 * point of view and with respect to the computer representation. Assignment to a larger trig_array
 * is allowed, assignment to smaller results in assertion failure.
 * @param[in] t2 right-hand side piranha::trig_array.
 */
      void assignment(const trig_array &t2)
      {
        increase_size(t2.g_width());
        private_container_=t2.private_container_;
      }
// Data members.
    private:
      container_type        private_container_;
  };

/// Resizable type-traits specialization for piranha::trig_array.
  template <>
    template <int Bits>
    struct is_resizable<trig_array<Bits> >:public true_type {};
}
#endif
