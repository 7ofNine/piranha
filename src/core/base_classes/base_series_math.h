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

namespace piranha
{
  // Do not use this to merge with self, assertion will fail.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool Sign, class Derived2, class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::merge_terms(const Derived2 &s2, const ArgsTuple &args_tuple)
  {
    p_assert(derived_cast != &s2);
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived2::const_sorted_iterator const_sorted_iterator2;
    const_sorted_iterator it_hint = derived_const_cast->template nth_index<0>().end();
    const const_sorted_iterator2 it_f = s2.template nth_index<0>().end();
    for (const_sorted_iterator2 it = s2.template nth_index<0>().begin(); it != it_f; ++it)
    {
      // No need to check, we are merging from another series.
      it_hint = insert<false,Sign>(*it,args_tuple,it_hint);
    }
  }
}

#endif
