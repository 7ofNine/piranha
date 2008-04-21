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

#include "../exceptions.h"

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
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::add_series(const Derived2 &s2)
  {
    return merge_with_series<true>(s2);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::subtract_series(const Derived2 &s2)
  {
    return merge_with_series<false>(s2);
  }

  // Multiply by generic entity. This means that the coefficients of the series get multiplied one by one
  // by the input value.
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class T>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::mult_by_generic(const T &x)
  {
    Derived retval;
    retval.merge_args(*derived_const_cast);
    derived_cast->multiply_coefficients_by(x,retval,m_arguments);
    swap(retval);
    return *derived_cast;
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class T>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::divide_by_generic(const T &x)
  {
    Derived retval;
    retval.merge_args(*derived_const_cast);
    derived_cast->divide_coefficients_by(x,retval,m_arguments);
    swap(retval);
    return *derived_cast;
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::mult_by_series(const Derived2 &s2)
  {
    Derived retval;
    // First we merge the arguments of the two series.
    merge_args(s2);
    // Then we assign the merged arguments sets to the return value series.
    retval.merge_args(*derived_const_cast);
    // Now we can multiply.
    derived_cast->multiply_by_series(s2,retval,m_arguments);
    // Swap terms with those accumulated into retval.
    derived_cast->swap_terms(retval);
    return *derived_cast;
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator+=(const int &n)
  {
    return derived_cast->template merge_with_number<true>(n,m_arguments);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator+=(const double &x)
  {
    return derived_cast->template merge_with_number<true>(x,m_arguments);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator+=(const Derived &s2)
  {
    return add_series(s2);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator-=(const int &n)
  {
    return derived_cast->template merge_with_number<false>(n,m_arguments);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator-=(const double &x)
  {
    return derived_cast->template merge_with_number<false>(x,m_arguments);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator-=(const Derived &s2)
  {
    return subtract_series(s2);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator*=(const int &n)
  {
    return mult_by_generic(n);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator*=(const double &x)
  {
    return mult_by_generic(x);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator*=(const Derived &s2)
  {
    return mult_by_series(s2);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator/=(const int &n)
  {
    return divide_by_generic(n);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator/=(const double &x)
  {
    return divide_by_generic(x);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::pow(const double &x) const throw(unsuitable)
  {
    Derived retval(derived_const_cast->a_pow(x,derived_const_cast->m_arguments));
    retval.m_arguments = derived_const_cast->m_arguments;
    return retval;
  }
}

#endif
