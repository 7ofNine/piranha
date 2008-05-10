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

#ifndef PIRANHA_CELMEC_TOOLBOX_H
#define PIRANHA_CELMEC_TOOLBOX_H

#include <string>

#include "../integer_typedefs.h"
#include "../psym.h"

namespace piranha
{
  /// Toolbox for Celestial Mechanics.
  /**
   * This toolbox requires the following toolboxes: multiplication toolbox, Poisson series toolbox,
   * special functions toolbox. Derived must be a Poisson series with
   * polynomial arguments in slot 0 of the arguments tuple.
   */
  template <class Derived>
    class celmec_toolbox
  {
    public:
      static Derived r_a(const Derived &e_series, const Derived &M_series)
      {
        // First let's build 1+1/2 e^2.
        Derived retval(e_series);
        retval *= e_series;
        retval /= (max_fast_int)2;
        retval += (max_fast_int)1;
        // Let's find out the upper limit of the r_a development, according to the truncation limits
        // set in the truncator. Minimum-exponent-wise, r_a is equivalent to a power series from 0 to infty
        // of e^(1+i).
        const size_t n = Derived::multiplier_type::truncator_type::power_series_limit(e_series,e_series.m_arguments);
        Derived tmp;
        for (size_t i = 1; i <= n; ++i)
        {
          Derived expansion_term(e_series);
          expansion_term *= (max_fast_int)i;
          expansion_term = expansion_term.dbesselJ((max_fast_int)i);
          expansion_term /= (max_fast_int)i;
          Derived trig(M_series);
          trig *= (max_fast_int)i;
          trig = trig.cos();
          expansion_term *= trig;
          tmp += expansion_term;
        }
        tmp *= e_series;
        tmp *= (max_fast_int)(-2);
        retval += tmp;
        return retval;
      }
      static Derived r_a(const psym &e, const psym &M) {return r_a(Derived(e),Derived(M));}
      static Derived r_a(const std::string &e_name, const std::string &M_name)
      {
        return r_a(*psym_manager::get_pointer(e_name),*psym_manager::get_pointer(M_name));
      }
      static Derived sin_f(const Derived &e_series, const Derived &M_series)
      {
        Derived tmp(e_series);
        tmp *= e_series;
        tmp = (max_fast_int)1 - tmp;
        tmp = tmp.pow(.5);
        tmp *= (max_fast_int)2;
        Derived retval;
        const size_t n = Derived::multiplier_type::truncator_type::power_series_limit(e_series,e_series.m_arguments);
        // Here we reach n+1 because the exponent starts from power 0, while i starts from 1.
        for (size_t i = 1; i <= (n+1); ++i)
        {
          Derived expansion_term(e_series);
          expansion_term *= (max_fast_int)i;
          expansion_term = expansion_term.dbesselJ((max_fast_int)i);
          Derived trig(M_series);
          trig *= (max_fast_int)i;
          trig = trig.sin();
          expansion_term *= trig;
          retval += expansion_term;
        }
        retval *= tmp;
        return retval;
      }
      static Derived sin_f(const psym &e, const psym &M) {return sin_f(Derived(e),Derived(M));}
      static Derived sin_f(const std::string &e_name, const std::string &M_name)
      {
        return sin_f(*psym_manager::get_pointer(e_name),*psym_manager::get_pointer(M_name));
      }
  };
}

#endif
