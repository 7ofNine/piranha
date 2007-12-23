/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_BASE_PSERIES_GETTERS_AND_SETTERS_H
#define PIRANHA_BASE_PSERIES_GETTERS_AND_SETTERS_H

#include "base_pseries_ta_macros.h"

namespace piranha
{
// GETTERS
// ===============================

/// Return series length.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline size_t base_pseries<__PIRANHA_BASE_PS_TP>::length() const
  {
    return g_series_set()->size();
  }

/// Return trigonometric width.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline size_t base_pseries<__PIRANHA_BASE_PS_TP>::trig_width() const
  {
//p_assert(arguments().template get<1>().size()==lin_args_.size());
    return arguments().template get<1>().size();
  }

/// Return coefficient width.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline size_t base_pseries<__PIRANHA_BASE_PS_TP>::cf_width() const
  {
    return arguments().template get<0>().size();
  }

/// Return const reference to the vector of linear arguments.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const vector_int16 &base_pseries<__PIRANHA_BASE_PS_TP>::lin_args() const
  {
    return lin_args_;
  }

/// Return const reference to the vector of coefficient arguments.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const vector_psym_p &
    base_pseries<__PIRANHA_BASE_PS_TP>::cf_args() const
  {
    return m_arguments.template get<0>();
  }

/// Return const reference to the vector of trigonometric arguments.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const vector_psym_p &
    base_pseries<__PIRANHA_BASE_PS_TP>::trig_args() const
  {
    return m_arguments.template get<1>();
  }

/// Return const reference to the tuple of arguments vectors.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const typename base_pseries<__PIRANHA_BASE_PS_TP>::arguments_tuple_type &
    base_pseries<__PIRANHA_BASE_PS_TP>::arguments() const
  {
    return m_arguments;
  }

/// Return const reference to the set of terms.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const typename base_pseries<__PIRANHA_BASE_PS_TP>::series_set_type *
    base_pseries<__PIRANHA_BASE_PS_TP>::g_series_set() const
  {
    return &private_series_set_;
  }

/// Return a const reference to the sorted index.
/**
 * @see base_pseries::sorted_index.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const typename base_pseries<__PIRANHA_BASE_PS_TP>::sorted_index &
    base_pseries<__PIRANHA_BASE_PS_TP>::g_s_index() const
  {
    return g_series_set()->template get<0>();
  }

/// Return a const reference to the hashed index.
/**
 * @see base_pseries::hashed_index.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline const typename base_pseries<__PIRANHA_BASE_PS_TP>::hashed_index &
    base_pseries<__PIRANHA_BASE_PS_TP>::g_h_index() const
  {
    return g_series_set()->template get<1>();
  }

/// Return index of argument "name".
/**
 * N is the index of arguments vector in the arguments vectors tuple. Return value is a
 * std::pair<bool,size_t> in which the first element reports if the argument was found or not,
 * while the second elemnt reports the index of the element.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <int N>
    std::pair<bool,size_t> base_pseries<__PIRANHA_BASE_PS_TP>::arg_index(const std::string &name) const
  {
    std::pair<bool,size_t> retval(false,0);
    const size_t w=arguments().template get<N>().size();
    for (size_t i=0;i<w;++i)
    {
      if (arguments().template get<N>()[i]->name() == name)
      {
        retval.first = true;
        retval.second =i;
        break;
      }
    }
    return retval;
  }

/// Return index of coefficient argument "name".
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline std::pair<bool,size_t> base_pseries<__PIRANHA_BASE_PS_TP>::cf_arg_index(const std::string &name) const
  {
    return arg_index<0>(name);
  }

/// Return index of trigonometric argument "name".
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline std::pair<bool,size_t> base_pseries<__PIRANHA_BASE_PS_TP>::trig_arg_index(const std::string &name) const
  {
    return arg_index<1>(name);
  }

/// Return a numerical value corresponding to the memory address of the series.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline size_t base_pseries<__PIRANHA_BASE_PS_TP>::address() const
  {
    return (size_t)(void *)this;
  }

/// Begin of series.
/**
 * Returns an iterator pointing to the first term of the series. This allows to mimic the behaviour
 * of an STL container.
 * @see base_pseries::end().
 * @see base_pseries::iterator.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index
    base_pseries<__PIRANHA_BASE_PS_TP>::begin() const
  {
    return g_s_index().begin();
  }

/// End of series.
/**
 * Returns an iterator pointing to the last term of the series. This allows to mimic the behaviour
 * of an STL container.
 * @see base_pseries::begin().
 * @see base_pseries::iterator.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index
    base_pseries<__PIRANHA_BASE_PS_TP>::end() const
  {
    return g_s_index().end();
  }

// SETTERS
// ===============================

/// Return mutable reference to the vector of linear arguments.
 template <__PIRANHA_BASE_PS_TP_DECL>
    inline vector_int16 &base_pseries<__PIRANHA_BASE_PS_TP>::lin_args()
  {
    return lin_args_;
  }

/// Return mutable reference to the vector of coefficient arguments.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline vector_psym_p &
    base_pseries<__PIRANHA_BASE_PS_TP>::cf_args()
  {
    return m_arguments.template get<0>();
  }

/// Return mutabl reference to the vector of trigonometric arguments.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline vector_psym_p &
    base_pseries<__PIRANHA_BASE_PS_TP>::trig_args()
  {
    return m_arguments.template get<1>();
  }

/// Return mutable reference to the tuple of arguments vectors.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::arguments_tuple_type &
    base_pseries<__PIRANHA_BASE_PS_TP>::arguments()
  {
    return m_arguments;
  }

/// Return mutable reference to the set of terms.
 template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::series_set_type *
    base_pseries<__PIRANHA_BASE_PS_TP>::s_series_set()
  {
    return &private_series_set_;
  }

/// Return mutable reference to the sorted index.
 template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::sorted_index &
    base_pseries<__PIRANHA_BASE_PS_TP>::s_s_index()
  {
    return s_series_set()->template get<0>();
  }

/// Return mutable reference to the hashed index.
 template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::hashed_index &
    base_pseries<__PIRANHA_BASE_PS_TP>::s_h_index()
  {
    return s_series_set()->template get<1>();
  }
}
#endif
