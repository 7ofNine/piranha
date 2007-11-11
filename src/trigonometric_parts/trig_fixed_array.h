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

#include "base_trig_array.h"
#include "common_typedefs.h"
#include "packed_int_array.h"

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
    class trig_fixed_array: public base_trig_array<Bits,trig_fixed_array<Dim,Bits> >
  {
      typedef base_trig_array<Bits,trig_fixed_array<Dim,Bits> > ancestor;
      template <int Bits_, class Derived>
        friend class base_trig_array;
    public:
      typedef packed_int_array<Dim,Bits> container_type;
      typedef typename ancestor::value_type value_type;
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_fixed_array():ancestor::base_trig_array()
        {}
/// Ctor from piranha::deque_string.
      trig_fixed_array(const deque_string &sd):ancestor::base_trig_array()
      {
// TODO: check for integer limits when building.
        const size_t w=sd.size();
        if (w == 0)
        {
          std::cout << "Warning: not enough elements to construct trig_fixed_array." << std::endl;
          std::abort();
          return;
        }
        const trig_size_t d=g_width();
// Now w >= 1.
        if ((w-1) != d)
        {
          std::cout << "Warning: wrong size for trig_fixed_array in ctor from string." << std::endl;
// TODO: Here we continue really, just remains as debug.
          std::abort();
        }
        trig_size_t i;
        for (i=0;i<boost::minmax(d,(trig_size_t)(w-1)).get<0>();++i)
        {
// TODO: lexical conversion from int: store in temp and check for boundaries.
          private_container_[i]=utils::lexical_converter<int>(sd[i]);
        }
        for (;i<d;++i)
        {
          private_container_[i]=0;
        }
// Take care of flavour.
        if (*sd.back().c_str()=='s')
        {
          ancestor::s_flavour()=false;
        }
      }
      ~trig_fixed_array()
        {}
// Manip.
      void prepend_args(const size_t &)
      {
        p_assert(false);
      }
      void append_args(const size_t &)
      {
        p_assert(false);
      }
      void increase_size(const size_t &w)
      {
        p_assert(w == g_width());
      }
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
// FIXME: do a real check here, cf_width in series could be different than 0.
      bool checkup(const size_t &) const
      {
        return true;
      }
      bool operator==(const trig_fixed_array &t2) const
      {
        if (ancestor::g_flavour() != t2.g_flavour())
        {
          return false;
        }
        return (private_container_ == t2.private_container_);
      }
      bool operator<(const trig_fixed_array &t2) const
      {
        return ancestor::less_than(t2);
      }
      size_t hasher() const
      {
        size_t seed=ancestor::g_flavour();
        private_container_.hasher(seed);
        return seed;
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
        container_type::trigmult(private_container_,t2.private_container_,
          ret1.private_container_,ret2.private_container_);
#else
        tfa_unrollers<dimension>::mult(&private_container_[0]+dimension,&t2.private_container_[0]+dimension,
          &ret1.private_container_[0]+dimension,&ret2.private_container_[0]+dimension);
#endif
      }
      trig_fixed_array &operator*=(const int &n)
      {
        ancestor::mult_by_int(n);
        return *this;
      }
// End INTERFACE definition.
//-------------------------------------------------------
    private:
      static const trig_size_t &g_width()
      {
        return dimension;
      }
      const value_type *g_container() const
      {
        return &private_container_[0];
      }
      value_type *s_container()
      {
        return &private_container_[0];
      }
// Data members.
    private:
      container_type            private_container_;
      static const trig_size_t  dimension = (trig_size_t)(Dim);
  };

  template <int Dim, int Bits>
    const trig_size_t trig_fixed_array<Dim,Bits>::dimension;
}

#endif
