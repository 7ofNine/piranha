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
#include <complex>

#include "../exceptions.h"
#include "../mp.h"
#include "../null_type.h"
#include "base_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Return true if series is a singleton coefficient.
	/**
	 * @return true if series has length one and the key of the only term is equivalent to unity, false otherwise.
	 */
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
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::eval_type
		toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_eval(const double &t, const ArgsTuple &args_tuple) const
	{
		const const_iterator it_f = end();
		eval_type retval(0);
		for (const_iterator it = begin(); it != it_f; ++it) {
			eval_type tmp(it->m_cf.eval(t, args_tuple));
			tmp *= it->m_key.eval(t, args_tuple);
			retval += tmp;
		}
		return retval;
	}

	/// Length of the series.
	/**
	 * @return number of terms in the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::size_type toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::length() const
	{
		return m_container.size();
	}

	/// Is series empty?
	/**
	 * @return true if the series has zero terms, false otherwise.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::empty() const
	{
		return m_container.empty();
	}

	/// Number of atoms in the series.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::size_type toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::atoms() const
	{
		size_type retval = 0;
		const const_iterator it_f = end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			retval += it->m_cf.atoms() + it->m_key.atoms();
		}
		return retval;
	}

	/// Base series equality test.
	/**
	 * Please note that this method has no knowledge about arguments: all comparisons performed here on coefficients and keys
	 * assume that the arguments tuples of this and other have been merged.
	 *
	 * @param[in] other series this will be compared to.
	 *
	 * @return false if: lengths of series differ or at least one term of this series is not found in other, true otherwise.
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

	// Helper functor to check generically if a number, scalar or complex, is zero.
	template <class T>
	struct numerical_comparison_zero_check
	{
		bool operator()(const T &value) const
		{
			return value == 0;
		}
	};

	template <class T>
	struct numerical_comparison_zero_check<std::complex<T> >
	{
		bool operator()(const std::complex<T> &value) const
		{
			return value.real() == 0 && value.imag() == 0;
		}
	};

	/// Base generic numerical equality test.
	/**
	 * Will return true in the following cases:
	 *
	 * - input argument is zero and series is empty,
	 * - series is singleton coefficient equal to x.
	 *
	 * Otherwise, false will be returned.
	 *
	 * @param[in] x numerical quantity.
	 *
	 * @return true if series is mathematically equivalent to the input numerical quantity.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::generic_numerical_comparison(const Number &x) const
	{
		if (numerical_comparison_zero_check<Number>()(x)) {
			return empty();
		}
		if (!is_single_cf()) {
			return false;
		}
		return (begin()->m_cf == x);
	}

	/// Base numerical equality test.
	/**
	 * @param[in] x argument of comparison.
	 *
	 * @return generic_numerical_comparison() on x.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const double &x) const
	{
		return generic_numerical_comparison(x);
	}

	/// Base numerical equality test.
	/**
	 * @param[in] q argument of comparison.
	 *
	 * @return generic_numerical_comparison() on q.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const mp_rational &q) const
	{
		return generic_numerical_comparison(q);
	}

	/// Base numerical equality test.
	/**
	 * @param[in] z argument of comparison.
	 *
	 * @return generic_numerical_comparison() on z.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_equal_to(const mp_integer &z) const
	{
		return generic_numerical_comparison(z);
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
