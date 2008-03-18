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

#include <valarray>

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
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class T, class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::multiply_coefficients_by(const T &x, Derived &retval,
    const ArgsTuple &args_tuple) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived::term_type term_type;
    // Make sure we are inserting into an empty return value.
    p_assert(retval.template nth_index<0>().empty());
    const_sorted_iterator it_hint = retval.template nth_index<0>().end();
    const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
    for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it)
    {
      term_type term(*it);
      term.m_cf.mult_by(x,args_tuple);
      it_hint = retval.insert(term,args_tuple,it_hint);
    }
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::cache_terms_pointers(std::valarray<term_type *> &v) const
    {
      typedef typename Derived::const_sorted_iterator const_sorted_iterator;
      const size_t l = derived_const_cast->template nth_index<0>().size();
      v.resize(l);
      size_t i=0;
      const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();;
      for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it)
      {
        v[i]=&(*it);
        ++i;
      }
    }

  // Multiply term-by-term with another series, and place the result into retval.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class Derived2, class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::multiply_by_series(const Derived2 &s2, Derived &retval,
    const ArgsTuple &args_tuple) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived2::const_sorted_iterator const_sorted_iterator2;
    // Make sure we are inserting into an empty return value.
    p_assert(retval.template nth_index<0>().empty());
    // Just leave an empty series if this or s2 are zero.
    if (derived_const_cast->template nth_index<0>.empty() or s2.template nth_index<0>.empty())
    {
      return;
    }
    // Optimize if the second series is a pure coefficient series.
    // TODO: test the effectiveness of this by multiplying with single cf series in the first and second place.
    if (s2.is_single_cf())
    {
      multiply_coefficients_by(s2.template nth_index<0>().begin()->m_cf,retval,args_tuple);
      return;
    }
    // Let's cache the iterators to the terms of the series into two separate vectors, in order to
    // speed up further manipulations.
    std::valarray<term_type *> vpt1, vpt2;
    derived_const_cast->cache_terms_pointers(vpt1);
    s2.cache_terms_pointers(vpt2);
    // Create the truncator class.
    Truncator<Derived,Derived2> trunc(*derived_const_cast,s2);
  }
}

#endif
