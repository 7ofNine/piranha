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

#ifndef PIRANHA_BASE_SERIES_PROBE_H
#define PIRANHA_BASE_SERIES_PROBE_H

#include <cmath> // For std::abs.
#include <cstddef>

#include "../exceptions.h"
#include "../mp.h"
#include "../null_type.h"
#include "base_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::is_single_cf() const
	{
		return (length() == 1 && begin()->m_key.is_unity());
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline double toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_norm(const ArgsTuple &args_tuple) const
	{
		const const_iterator it_f = end();
		double retval = 0;
		for (const_iterator it = begin(); it != it_f; ++it) {
			retval += it->m_cf.norm(args_tuple) * it->m_key.norm(args_tuple);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_eval_type
		toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_eval(const double &t, const ArgsTuple &args_tuple) const
	{
		const const_iterator it_f = end();
		base_eval_type retval(0);
		for (const_iterator it = begin(); it != it_f; ++it) {
			base_eval_type tmp(it->m_cf.eval(t, args_tuple));
			tmp *= it->m_key.eval(t, args_tuple);
			retval += tmp;
		}
		return retval;
	}

	/// Return the number of terms of the series.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline std::size_t toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::length() const
	{
		return m_container.size();
	}

	/// Is series empty?
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::empty() const
	{
		return m_container.empty();
	}

	/// Number of atoms in the series.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline std::size_t toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::atoms() const
	{
		std::size_t retval = 0;
		const const_iterator it_f = end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			retval += it->m_cf.atoms() + it->m_key.atoms();
		}
		return retval;
	}

	/// Test for equality.
	/**
	 * Please note that this method has no knowledge about arguments: all comparisons performed here on coefficients and keys
	 * assume that the arguments tuples of this and other have been merged.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const Derived &other) const
	{
		if (length() != other.length()) {
			return false;
		}
		const const_iterator it_f = end(), it_f_other = other.end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			const_iterator it_other(other.find_term(*it));
			if (it_other == it_f_other || !(it_other->m_cf == it->m_cf)) {
				return false;
			}
		}
		return true;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::generic_numerical_comparison(const Number &x) const
	{
		// Use std::abs() to cope with complex numbers. In case of mp classes this could throw value_error,
		// so handle this.
		try {
			if (std::abs(x) == 0) {
				return empty();
			}
		} catch (const value_error &) {
			// Don't do anything, continue the flow.
		}
		if (!is_single_cf()) {
			return false;
		}
		return (begin()->m_cf == x);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const double &x) const
	{
		return generic_numerical_comparison(x);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const mp_rational &q) const
	{
		return generic_numerical_comparison(q);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const mp_integer &z) const
	{
		return generic_numerical_comparison(z);
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
