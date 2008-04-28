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

#include "../base_classes/truncators.h" // To establish a limit for the binomial expansion.
#include "../exceptions.h"
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
        Derived retval;
        const int pow_n = (int)nearbyint(y);
        // If y is an integer, assert that it is -1 (i.e., that we are inverting).
        const bool is_inversion = std::abs(y - pow_n) < settings::numerical_zero();
        p_assert(!is_inversion or pow_n == -1);
        term_type A(*derived_const_cast->template nth_index<0>().begin());
        // If we are dealing with a single monomial and the power is -1, just try to invert it.
        if (is_inversion and derived_const_cast->template nth_index<0>().size() == 1)
        {
          term_type tmp(A);
          tmp.m_cf = tmp.m_cf.pow(-1,args_tuple);
          tmp.m_key.invert_sign();
          retval.insert(tmp,args_tuple,retval.template nth_index<0>().end());
        }
        else
        {
          // Start the binomial expansion.
          term_type tmp_term;
          // First we calculate A**y. See if we can raise to real power the coefficient.
          tmp_term.m_cf = A.m_cf.pow(y,args_tuple);
          // If we are inverting, it is all fine and dandy.
          if (is_inversion)
          {
            tmp_term.m_key = A.m_key;
            tmp_term.m_key.invert_sign();
          }
          // Otherwise the only case in which we can proceed is when the key is equal to 1.
          else if (!A.m_key.is_unity())
          {
            // Can't handle real exponentiation of non-unitary sets of exponents.
            throw (unsuitable("Cannot raise to real power leading monomial during binomial expansion of polynomial."));
          }
          Derived Apowy;
          Apowy.insert(tmp_term,args_tuple,Apowy.template nth_index<0>().end());
          // Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
          tmp_term.m_cf = A.m_cf.pow(-1,args_tuple);
          tmp_term.m_key = A.m_key;
          tmp_term.m_key.invert_sign();
          Derived Ainv;
          Ainv.insert(tmp_term,args_tuple,Ainv.template nth_index<0>().end());
          // Now, on to X/A.
          // This is X, i.e., the original series without the leading term, which will then be divided by A.
          Derived XoverA(*derived_const_cast);
          XoverA.template term_erase<0>(args_tuple,XoverA.template nth_index<0>().begin());
          XoverA.mult_by(Ainv,args_tuple);
          if (XoverA.empty())
          {
            throw unsuitable("The argument of the binomial expansion is empty.");
          }
          const size_t n(binomial_limit(XoverA,args_tuple));
          __PDEBUG(std::cout << "Calculated limit for binomial expansion: " << n << '\n');
          // Let's proceed now to the bulk of the binomial expansion. Luckily we can compute the needed generalised
          // binomial coefficient incrementally at every step. We start with 1.
          Derived tmp(1,args_tuple);
          retval.add(tmp,args_tuple);
          for (size_t i = 1; i <= n; ++i)
          {
            tmp.mult_by(y-i+1,args_tuple);
            tmp.divide_by((int)i,args_tuple);
            tmp.mult_by(XoverA,args_tuple);
            retval.add(tmp,args_tuple);
          }
          // Finally, multiply the result of the summation by A**y.
          retval.mult_by(Apowy,args_tuple);
        }
        return retval;
      }
    private:
      template <class ArgsTuple>
        static size_t binomial_limit(const Derived &XoverA, const ArgsTuple &args_tuple)
      {
        typedef typename Derived::term_type term_type;
        typedef typename Derived::const_sorted_iterator const_sorted_iterator;
        // First let's translate the psym_p - integer pairs from the exponent truncator into
        // size_t - integer pairs. I.e., establish which truncation limits are relevant
        // to the given tuple of arguments sets and thei positions.
        const std::vector<std::pair<size_t,int> > pos(base_expo_truncator::get_positions(
          args_tuple.template get<term_type::key_type::position>()));
        const size_t pos_size = pos.size();
        if (pos_size == 0)
        {
          throw not_existing("Cannot establish the limit of the binomial expansion: the limits set in the "
           "exponent truncator do not involve arguments belonging to this series.");
        }
        // We don't want negative exponent limits, because they lead to a never-ending expansion.
        for (size_t i = 0; i < pos_size; ++i)
        {
          if (pos[i].second < 0)
          {
            throw unsuitable("Cannot proceed to binomial expansion if there are negative exponent limits for this series.");
          }
        }
        // The code that follows is used to build a vector of minimum degrees for those symbolic
        // arguments subject to exponent truncation.
        std::vector<int> minimum_degrees(pos_size);
        // We must do this the first time outside of the for cycle to set the initial values.
        const_sorted_iterator it = XoverA.template nth_index<0>().begin();
        for (size_t i = 0; i < pos_size; ++i)
        {
          minimum_degrees[i] = it->m_key.degree_of(pos[i].first);
        }
        ++it;
        int tmp;
        const const_sorted_iterator it_f = XoverA.template nth_index<0>().end();
        for (; it != it_f; ++it)
        {
          for (size_t i = 0; i < pos_size; ++i)
          {
            tmp = it->m_key.degree_of(pos[i].first);
            if (tmp < minimum_degrees[i])
            {
              minimum_degrees[i] = tmp;
            }
          }
        }
        // Calculate the binomial expansion limit given the minimum degrees we just calculated.
        std::vector<size_t> bin_limits(pos_size);
        for (size_t i = 0; i < pos_size; ++i)
        {
          // If the minimum degree of the symbol whose exponent we want to limit is zero or less, we are
          // screwed, since the degree of the symbol would not increase at every step of the expansion
          // and we would end up in an infinite loop.
          if (minimum_degrees[i] <= 0)
          {
            throw unsuitable("An infinite binomial expansion would be needed to satisfy the given limits on exponents.");
          }
          bin_limits[i] =(size_t)std::ceil((float)pos[i].second / minimum_degrees[i]);
        }
        // Finally, the biggest limit is the one that defines the length of the binomial expansion.
        return *(std::max_element(bin_limits.begin(),bin_limits.end()));
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
