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

#ifndef PIRANHA_TRIG_SIMD_ARRAY_H
#define PIRANHA_TRIG_SIMD_ARRAY_H

#include <boost/algorithm/minmax.hpp>

#include "base_trig_array.h"
#include "simd_array.h"

namespace piranha
{
/// Trigonometric array, fixed size SSE2-enabled version.
  template <int Dim>
    class trig_simd_array: public base_trig_array<trig_simd_array<Dim> >
  {
      typedef base_trig_array<trig_simd_array<Dim> > ancestor;
      template <class Derived>
        friend class base_trig_array;
      typedef simd_array<Dim> container_type;
    public:
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_simd_array():ancestor::base_trig_array(),private_container_()
        {}
/// Copy ctor.
      trig_simd_array(const trig_simd_array &t):ancestor::base_trig_array(t),private_container_(t.private_container_)
        {}
/// Ctor from piranha::deque_string.
      trig_simd_array(const deque_string &sd):ancestor::base_trig_array()
      {
        const size_t w=sd.size();
        if (w == 0)
        {
          std::cout << "Warning: not enough elements to construct trig_simd_array." << std::endl;
          std::abort();
          return;
        }
        const usint d=g_width();
// Now w >= 1.
        if ((w-1) != d)
        {
          std::cout << "Warning: wrong size for trig_simd_array in ctor from string." << std::endl;
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
      ~trig_simd_array()
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
      void invert_sign()
      {
        private_container_.invert_sign();
      }
// Probing.
/// Data footprint.
/**
 * Returns the memory occupied by the data members.
 */
      size_t data_footprint() const
      {
        return sizeof(trig_simd_array);
      }
      bool checkup(const size_t &) const
      {
        return true;
      }
      bool operator==(const trig_simd_array &t2) const
      {
        if (ancestor::g_flavour() != t2.g_flavour())
        {
          return false;
        }
        return private_container_.equality_test(t2.private_container_);
      }
      bool operator<(const trig_simd_array &t2) const
      {
        if (ancestor::g_flavour() < t2.g_flavour())
        {
          return true;
        }
        else if (ancestor::g_flavour() > t2.g_flavour())
        {
          return false;
        }
        return private_container_.less_than(t2.private_container_);
      }
      size_t hasher() const
      {
        size_t seed=ancestor::g_flavour();
        private_container_.hasher(seed);
        return seed;
      }
      short int sign() const
      {
        return private_container_.sign();
      }
      bool is_zero() const
      {
        return private_container_.is_zero();
      }
// Math.
/// Multiplication.
      void trigmult(const trig_simd_array &t2, trig_simd_array &ret1, trig_simd_array &ret2) const
      {
        ret1.private_container_=private_container_;
        ret1.private_container_-=t2.private_container_;
        ret2.private_container_=private_container_;
        ret2.private_container_+=t2.private_container_;
      }
      trig_simd_array &operator=(const trig_simd_array &t2)
      {
        ancestor::assignment_operator(t2);
        return *this;
      }
      trig_simd_array &operator*=(const mult_t &n)
      {
        private_container_*=n;
        return *this;
      }
// End INTERFACE definition.
//-------------------------------------------------------
    private:
      static const usint &g_width()
      {
        return dimension;
      }
      const container_type &g_container() const
      {
        return private_container_;
      }
      container_type &s_container()
      {
        return private_container_;
      }
      void assignment(const trig_simd_array &t2)
      {
        private_container_=t2.private_container_;
      }
// Data members.
    private:
      container_type      private_container_;
      static const usint  dimension = (usint)Dim;
  };

  template <int Dim>
    const usint trig_simd_array<Dim>::dimension;
}

#endif
