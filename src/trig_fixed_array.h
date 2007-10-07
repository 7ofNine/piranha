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

#include <boost/static_assert.hpp>
#include <boost/algorithm/minmax.hpp>
#include <cstring>

#include "base_trig_array.h"

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

  template <int N>
    inline void mult_unroller(const int16 *end_array1, const int16 *end_array2, int16 *ret_array1, int16 *ret_array2)
  {
    ret_array1[-N]=end_array1[-N]-end_array2[-N];
    ret_array2[-N]=end_array1[-N]+end_array2[-N];
    mult_unroller<N-1>(end_array1,end_array2,ret_array1,ret_array2);
  }

  template <>
    inline void mult_unroller<1>(const int16 *end_array1, const int16 *end_array2, int16 *ret_array1, int16 *ret_array2)
  {
    ret_array1[-1]=end_array1[-1]-end_array2[-1];
    ret_array2[-1]=end_array1[-1]+end_array2[-1];
  }

// TODO: Are these two used?
//   template <int N>
//     inline void invert_unroller(int16 *end_array)
//   {
//     end_array[-N]=-end_array[-N];
//     invert_unroller<N-1>(end_array);
//   }
//
//   template <>
//     inline void invert_unroller<1>(int16 *end_array)
//   {
//     end_array[-1]=-end_array[-1];
//   }

/// Trigonometric array, fixed size version.
  template <int Dim>
    class trig_fixed_array: public base_trig_array<trig_fixed_array<Dim> >
  {
// Check that dimension is sane.
      BOOST_STATIC_ASSERT(Dim > 0);
      BOOST_STATIC_ASSERT(Dim < 100);
      typedef base_trig_array<trig_fixed_array<Dim> > ancestor;
      template <class Derived>
        friend class base_trig_array;
    public:
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_fixed_array():ancestor::base_trig_array()
        {}
/// Copy ctor.
      trig_fixed_array(const trig_fixed_array &t):ancestor::base_trig_array(t)
        {
          assignment(t);
        }
/// Ctor from piranha::deque_string.
      trig_fixed_array(const deque_string &sd):ancestor::base_trig_array()
      {
        const size_t w=sd.size();
        if (w == 0)
        {
          std::cout << "Warning: not enough elements to construct trig_fixed_array." << std::endl;
          std::abort();
          return;
        }
        const uint16 d=g_width();
// Now w >= 1.
        if ((w-1) != d)
        {
          std::cout << "Warning: wrong size for trig_fixed_array in ctor from string." << std::endl;
// TODO: Here we continue really, just remains as debug.
          std::abort();
        }
        size_t i;
        for (i=0;i<boost::minmax((size_t)d,w-1).get<0>();++i)
        {
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
        return ancestor::equality_test(t2);
      }
      bool operator<(const trig_fixed_array &t2) const
      {
        return ancestor::less_than(t2);
      }
      size_t hasher() const
      {
        size_t seed=ancestor::g_flavour();
        hash_unroller<dimension>(seed,private_container_+dimension);
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
        mult_unroller<dimension>(private_container_+dimension,t2.private_container_+dimension,
          ret1.private_container_+dimension,ret2.private_container_+dimension);
      }
      trig_fixed_array &operator=(const trig_fixed_array &t2)
      {
        ancestor::assignment_operator(t2);
        return *this;
      }
      trig_fixed_array &operator*=(const int16 &n)
      {
        ancestor::mult_by_int16(n);
        return *this;
      }
// End INTERFACE definition.
//-------------------------------------------------------
    private:
      static const uint16 &g_width()
      {
        return dimension;
      }
      const int16 *g_container() const
      {
        return private_container_;
      }
      int16 *s_container()
      {
        return private_container_;
      }
      void assignment(const trig_fixed_array &t2)
      {
        memcpy((void *)private_container_,(const void *)t2.private_container_,sizeof(int16)*g_width());
      }
// Data members.
    private:
      int16              private_container_[Dim];
      static const uint16  dimension = (uint16)Dim;
  };

  template <int Dim>
    const uint16 trig_fixed_array<Dim>::dimension;
}

#endif
