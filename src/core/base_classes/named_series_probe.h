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

#ifndef PIRANHA_NAMED_SERIES_PROBE_H
#define PIRANHA_NAMED_SERIES_PROBE_H

#include "../config.h" // For (un)likely

namespace piranha
{
  // TMP for arguments compatibility check.
  template <class ArgsTuple>
    inline bool named_series_is_args_compatible(const ArgsTuple &a1, const ArgsTuple &a2)
  {
    const size_t w = a2.get_head().size();
    if (unlikely(w > a1.get_head().size()))
    {
      return false;
    }
    for (size_t i=0; i<w; ++i)
    {
      if (unlikely(a1.get_head()[i] != a2.get_head()[i]))
      {
        return false;
      }
    }
    return named_series_is_args_compatible(a1.get_tail(),a2.get_tail());
  }

  template <>
    inline bool named_series_is_args_compatible<boost::tuples::null_type>(
    const boost::tuples::null_type &, const boost::tuples::null_type &)
  {
    return true;
  }

  /// Compatibility check for arguments.
  /**
   * Test whether series' arguments are compatible with those from ps2. Compatibility
   * means that the number of arguments in all arguments sets are equal to or greater than ps2's, and
   * that arguments have the same positions as in ps2's.
   * @param[in] ps2 series compatibility is tested against.
   */
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline bool named_series<__PIRANHA_NAMED_SERIES_TP>::is_args_compatible(const Derived2 &ps2) const
  {
    return named_series_is_args_compatible(m_arguments,ps2.m_arguments);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline double named_series<__PIRANHA_NAMED_SERIES_TP>::norm() const
  {
    return derived_const_cast->b_norm(m_arguments);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline typename named_series<__PIRANHA_NAMED_SERIES_TP>::eval_type
    named_series<__PIRANHA_NAMED_SERIES_TP>::eval(const double &t) const
  {
    return derived_const_cast->b_eval(t,m_arguments);
  }
}

#endif
