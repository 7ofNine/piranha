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

#ifndef PIRANHA_COMMON_POLYNOMIAL_TOOLBOX_H
#define PIRANHA_COMMON_POLYNOMIAL_TOOLBOX_H

#include <cmath> // For nearbyint.
#include <string>

#include "../base_classes/series_math.h" // For the binomial expansion.
#include "../exceptions.h"
#include "../p_assert.h"
#include "../settings.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  template <class Derived>
    struct common_polynomial_toolbox
  {
    protected:
      /// Real power.
      /**
      * This method is written to work in conjunction with base_series::b_pow.
      */
      template <class ArgsTuple>
        Derived real_pow(const double &y, const ArgsTuple &args_tuple) const
      {
        typedef typename Derived::term_type term_type;
        // Here we know that the cases of single coefficient, empty series and natural power have already been taken care of
        // in base_series::b_pow. We also know that if y is an integer, it must be -1.
        p_assert(!derived_const_cast->is_single_cf() and !derived_const_cast->empty());
        const int pow_n = (int)nearbyint(y);
        // If y is an integer, assert that it is -1 (i.e., that we are inverting).
        const bool is_inversion = std::abs(y - pow_n) < settings::numerical_zero();
        p_assert(!is_inversion or pow_n == -1);
        term_type A(*derived_const_cast->template nth_index<0>().begin());
        // If we are dealing with a single monomial and the power is -1, just try to invert it.
        if (is_inversion and derived_const_cast->template nth_index<0>().size() == 1)
        {
          Derived retval;
          term_type tmp(A);
          tmp.m_cf = tmp.m_cf.pow(-1,args_tuple);
          tmp.m_key.invert_sign();
          retval.insert(tmp,args_tuple,retval.template nth_index<0>().end());
          return retval;
        }
        else
        {
          // This is X, i.e., the original series without the leading term, which will then be divided by A.
          Derived XoverA(*derived_const_cast);
          XoverA.template term_erase<0>(args_tuple,XoverA.template nth_index<0>().begin());
          // Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
          term_type tmp_term;
          tmp_term.m_cf = A.m_cf.pow(-1,args_tuple);
          tmp_term.m_key = A.m_key.pow(-1,args_tuple);
          Derived Ainv;
          Ainv.insert(tmp_term,args_tuple,Ainv.template nth_index<0>().end());
          // Now let's compute X/A.
          XoverA.mult_by(Ainv,args_tuple);
          // Get the expansion limit from the truncator.
          size_t n;
          try
          {
            n = Derived::multiplier_type::truncator_type::power_series_limit(XoverA,args_tuple);
          }
          catch (const unsuitable &u)
          {
            throw unsuitable(std::string("Polynomial is unsuitable for real exponentiation.\nThe reported error is: ")
              + u.what());
          }
          return binomial_expansion(A,XoverA,y,n,args_tuple);
        }
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
