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

#include <algorithm> // To calculate min degree.
#include <cmath>

// TODO: drop some of these includes later.
#include "../base_classes/expo_truncator.h" // To establish a limit for the binomial expansion.
#include "../base_classes/series_math.h" // For the binomial expansion.
#include "../exceptions.h"
#include "../p_assert.h"
#include "../settings.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  /// Common polynomial toolbox.
  /**
  * This toolbox assumes that the monomials are sorted in ascending total degree.
  */
  template <class Derived>
    struct common_polynomial_toolbox
  {
      /// Get the degree of the polynomial.
      int degree() const
      {
        if (derived_const_cast->template nth_index<0>().empty())
        {
          return 0;
        }
        typename Derived::const_sorted_iterator it = derived_const_cast->template nth_index<0>().end();
        --it;
        return it->m_key.degree();
      }
      /// Get the minimum degree of the polynomial.
      int min_degree() const
      {
        if (derived_const_cast->template nth_index<0>().empty())
        {
          return 0;
        }
        return derived_const_cast->template nth_index<0>().begin()->m_key.degree();
      }
//     protected:
//       /// Real power.
//       /**
//       * This method is written to work in conjunction with base_series::b_pow.
//       */
//       template <class ArgsTuple>
//         Derived real_pow(const double &y, const ArgsTuple &args_tuple) const
//       {
//         typedef typename Derived::term_type term_type;
//         // Here we know that the cases of single coefficient, empty series and natural power have already been taken care of
//         // in base_series::b_pow. We also know that if y is an integer, it must be -1.
//         p_assert(!derived_const_cast->is_single_cf() and !derived_const_cast->empty());
//         const int pow_n = (int)nearbyint(y);
//         // If y is an integer, assert that it is -1 (i.e., that we are inverting).
//         const bool is_inversion = std::abs(y - pow_n) < settings::numerical_zero();
//         p_assert(!is_inversion or pow_n == -1);
//         term_type A(*derived_const_cast->template nth_index<0>().begin());
//         // If we are dealing with a single monomial and the power is -1, just try to invert it.
//         if (is_inversion and derived_const_cast->template nth_index<0>().size() == 1)
//         {
//           Derived retval;
//           term_type tmp(A);
//           tmp.m_cf = tmp.m_cf.pow(-1,args_tuple);
//           tmp.m_key.invert_sign();
//           retval.insert(tmp,args_tuple,retval.template nth_index<0>().end());
//           return retval;
//         }
//         else
//         {
//           // This is X, i.e., the original series without the leading term, which will then be divided by A.
//           Derived XoverA(*derived_const_cast);
//           XoverA.template term_erase<0>(args_tuple,XoverA.template nth_index<0>().begin());
//           // Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
//           term_type tmp_term;
//           tmp_term.m_cf = A.m_cf.pow(-1,args_tuple);
//           tmp_term.m_key = A.m_key.pow(-1,args_tuple);
//           Derived Ainv;
//           Ainv.insert(tmp_term,args_tuple,Ainv.template nth_index<0>().end());
//           // Now let's compute X/A.
//           XoverA.mult_by(Ainv,args_tuple);
//           const size_t n = power_series_limit(XoverA,args_tuple);
//           __PDEBUG(std::cout << "Calculated limit for binomial expansion: " << n << '\n');
//           return binomial_expansion(A,XoverA,y,n,args_tuple);
//         }
//       }
//     private:
//       template <class ArgsTuple>
//         std::vector<int> min_limited_exponents(const ArgsTuple &args_tuple) const
//       {
//         std::vector<int> min_expos(derived_const_cast->min_key_ints());
        


//         typedef typename Derived::term_type term_type;
//         typedef typename Derived::const_sorted_iterator const_sorted_iterator;
//         // First let's translate the psym_p - integer pairs from the exponent truncator into
//         // size_t - integer pairs. I.e., establish which truncation limits are relevant
//         // to the given tuple of arguments sets and their positions.
//         const std::vector<std::pair<size_t,int> > pos(base_expo_truncator::get_positions_limits(
//           args_tuple.template get<term_type::key_type::position>()));
//         const size_t pos_size = pos.size();
//         std::vector<int> retval(pos_size);
//         // Just take a shortcut in case there are no relevant limited exponents.
//         if (pos_size == 0)
//         {
//           return retval;
//         }
//         // We don't want this to happen.
//         if (derived_const_cast->empty())
//         {
//           throw unsuitable("Cannot determine minimum limited exponents because series is empty.");
//         }
// //         // We don't want negative exponent limits, because they lead to a never-ending expansion.
// //         for (size_t i = 0; i < pos_size; ++i)
// //         {
// //           if (pos[i].second < 0)
// //           {
// //             throw unsuitable("Cannot proceed to a power series expansion if there are negative exponent limits.");
// //           }
// //         }
//         // We must do this the first time outside of the for cycle to set the initial values.
//         const_sorted_iterator it = XoverA.template nth_index<0>().begin();
//         for (size_t i = 0; i < pos_size; ++i)
//         {
//           retval[i] = it->m_key.degree_of(pos[i].first);
//         }
//         ++it;
//         int tmp;
//         const const_sorted_iterator it_f = XoverA.template nth_index<0>().end();
//         for (; it != it_f; ++it)
//         {
//           for (size_t i = 0; i < pos_size; ++i)
//           {
//             tmp = it->m_key.degree_of(pos[i].first);
//             if (tmp < retval[i])
//             {
//               retval[i] = tmp;
//             }
//           }
//         }
//         return retval;
//       }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
