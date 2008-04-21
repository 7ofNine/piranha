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

#ifndef PIRANHA_BASE_SERIES_MATH_H
#define PIRANHA_BASE_SERIES_MATH_H

#include <cmath>

#include "../exceptions.h"
#include "../settings.h"

namespace piranha
{
  // Do not use this to merge with self, assertion will fail.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool Sign, class Derived2, class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::merge_terms(const Derived2 &s2, const ArgsTuple &args_tuple)
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived2::const_sorted_iterator const_sorted_iterator2;
    p_assert(derived_cast != &s2);
    const_sorted_iterator it_hint = derived_const_cast->template nth_index<0>().end();
    const const_sorted_iterator2 it_f = s2.template nth_index<0>().end();
    for (const_sorted_iterator2 it = s2.template nth_index<0>().begin(); it != it_f; ++it)
    {
      // No need to check, we are merging from another series.
      it_hint = insert<false,Sign>(*it,args_tuple,it_hint);
    }
  }

  // Multiply all the coefficients of the series by a generic quantity x, and place the result into retval.
  // In case the coefficient is another series, the corresponding arguments must have been merged previously,
  // otherwise an assertion will fail when inserting terms.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class T, class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::multiply_coefficients_by(const T &x, Derived &retval,
    const ArgsTuple &args_tuple) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived::sorted_iterator sorted_iterator;
    typedef typename Derived::term_type term_type;
    // Make sure we are inserting into an empty return value.
    p_assert(retval.template nth_index<0>().empty());
    sorted_iterator it_hint = retval.template nth_index<0>().end();
    const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
    for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it)
    {
      term_type term(*it);
      term.m_cf.mult_by(x,args_tuple);
      it_hint = retval.insert(term,args_tuple,it_hint);
    }
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class T, class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::divide_coefficients_by(const T &x, Derived &retval,
    const ArgsTuple &args_tuple) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived::sorted_iterator sorted_iterator;
    typedef typename Derived::term_type term_type;
    // Make sure we are inserting into an empty return value.
    p_assert(retval.template nth_index<0>().empty());
    sorted_iterator it_hint = retval.template nth_index<0>().end();
    const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
    for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it)
    {
      term_type term(*it);
      term.m_cf.divide_by(x,args_tuple);
      it_hint = retval.insert(term,args_tuple,it_hint);
    }
  }

  /// Merge series with a number.
  /**
   * Term is constructed from coefficient constructed from number and default key, and then inserted into the series.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool Sign, class Number, class ArgsTuple>
    inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::merge_with_number(const Number &n, const ArgsTuple &args_tuple)
  {
    typename Derived::term_type term(typename Derived::term_type::cf_type(n,args_tuple),typename Derived::term_type::key_type());
    insert<true,Sign>(term,args_tuple,derived_cast->template nth_index<0>().end());
    return *derived_cast;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::mult_by(const int &n, const ArgsTuple &args_tuple)
  {
    Derived retval;
    multiply_coefficients_by(n,retval,args_tuple);
    swap_terms(retval);
    return *derived_cast;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::mult_by(const double &x, const ArgsTuple &args_tuple)
  {
    Derived retval;
    multiply_coefficients_by(x,retval,args_tuple);
    swap_terms(retval);
    return *derived_cast;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::mult_by(const Derived &s2, const ArgsTuple &args_tuple)
  {
    Derived retval;
    derived_cast->multiply_by_series(s2,retval,args_tuple);
    // Grab the terms accumulated into return value.
    swap_terms(retval);
    return *derived_cast;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::a_pow(const double &x,
    const ArgsTuple &args_tuple) const throw(unsuitable)
  {
    const int n = (int)nearbyint(x);
    if (std::abs(x - n) <= settings::numerical_zero())
    {
      if (n < 0)
      {
        Derived tmp(derived_const_cast->real_pow(-1,args_tuple));
        Derived retval(tmp.natural_pow((size_t)(-n),args_tuple));
        return retval;
      }
      else
      {
        Derived retval(natural_pow((size_t)n,args_tuple));
        return retval;
      }
    }
    else
    {
      Derived retval(derived_const_cast->real_pow(x,args_tuple));
      return retval;
    }
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::real_pow(const double &, const ArgsTuple &) const
  {
    throw (unsuitable("Unable to raise to real power."));
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::natural_pow(const size_t &n, const ArgsTuple &args_tuple) const
  {
    typedef typename Derived::term_type term_type;
    typedef typename term_type::cf_type cf_type;
    typedef typename term_type::key_type key_type;
    Derived retval;
    switch (n)
    {
      case 0:
      {
        retval.insert(term_type(cf_type(1,args_tuple),key_type()),args_tuple,derived_const_cast->template nth_index<0>().end());
        break;
      }
      case 1:
      {
        retval.m_container = derived_const_cast->m_container;
        break;
      }
      case 2:
      {
        retval.m_container = derived_const_cast->m_container;
        retval.mult_by(*derived_const_cast,args_tuple);
        break;
      }
      case 3:
      {
        retval.m_container = derived_const_cast->m_container;
        retval.mult_by(*derived_const_cast,args_tuple);
        retval.mult_by(*derived_const_cast,args_tuple);
        break;
      }
      case 4:
      {
        retval.m_container = derived_const_cast->m_container;
        retval.mult_by(*derived_const_cast,args_tuple);
        retval.mult_by(*derived_const_cast,args_tuple);
        retval.mult_by(*derived_const_cast,args_tuple);
        break;
      }
      // Exponentiation by squaring.
      default:
      {
        retval.insert(term_type(cf_type(1,args_tuple),key_type()),args_tuple,derived_const_cast->template nth_index<0>().end());
        Derived tmp(*derived_const_cast);
        size_t i = n;
        while (i)
        {
          if (i & 1)
          {
            retval.mult_by(tmp,args_tuple);
            --i;
          }
          tmp.mult_by(tmp,args_tuple);
          i/=2;
        }
      }
    }
    return retval;
  }
}

#endif
