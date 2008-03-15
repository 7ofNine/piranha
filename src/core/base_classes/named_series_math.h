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

#ifndef PIRANHA_NAMED_SERIES_MATH_H
#define PIRANHA_NAMED_SERIES_MATH_H

namespace piranha
{
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <bool Sign, class Derived2>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::merge_with_series(const Derived2 &s2)
  {
    // If we are merging with self, create a copy and call recursively.
    if ((void *)derived_cast == (void *)(&s2))
    {
      std::cout << "Merging with self, performing a copy." << std::endl;
      return merge_with_series<Sign,Derived2>(Derived(*derived_const_cast));
    }
    else
    {
      merge_args(s2);
      derived_cast->template merge_terms<Sign>(s2,m_arguments);
    }
    return *derived_cast;
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::add(const Derived2 &s2)
  {
    return merge_with_series<true>(s2);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::subtract(const Derived2 &s2)
  {
    return merge_with_series<false>(s2);
  }
}

#endif
