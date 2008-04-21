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

#ifndef PIRANHA_CF_SERIES_MATH_H
#define PIRANHA_CF_SERIES_MATH_H

namespace piranha
{
  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::add(const Derived &s2, const ArgsTuple &args_tuple)
  {
    derived_cast->template merge_terms<true>(s2,args_tuple);
    return *derived_cast;
  }

  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::subtract(const Derived &s2, const ArgsTuple &args_tuple)
  {
    derived_cast->template merge_terms<false>(s2,args_tuple);
    return *derived_cast;
  }

  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline void cf_series<__PIRANHA_CF_SERIES_TP>::invert_sign(const ArgsTuple &args_tuple)
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived::term_type term_type;
    Derived retval;
    const_sorted_iterator it_hint = retval.template nth_index<0>().end();
    const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
    for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it)
    {
      term_type term(*it);
      term.m_cf.invert_sign(args_tuple);
      // No need to check, we are merging terms from this series.
      it_hint = retval.template insert<false,true>(term,args_tuple,it_hint);
    }
    derived_cast->swap_terms(retval);
  }

  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class T, class ArgsTuple>
    inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divide_by_generic(const T &x, const ArgsTuple &args_tuple)
  {
    Derived retval;
    derived_cast->divide_coefficients_by(x,retval,args_tuple);
    derived_cast->swap_terms(retval);
    return *derived_cast;
  }

  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divide_by(const int &n, const ArgsTuple &args_tuple)
  {
    return divide_by_generic(n,args_tuple);
  }

  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divide_by(const double &x, const ArgsTuple &args_tuple)
  {
    return divide_by_generic(x,args_tuple);
  }
}

#endif