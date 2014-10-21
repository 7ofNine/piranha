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
#include "base_series_def.h"
#include "base_series_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Return true if series is a singleton coefficient.
	/**
	 * @return true if series has length one and the key of the only term is equivalent to unity, false otherwise.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool BaseSeries<__PIRANHA_BASE_SERIES_TP>::isSingleCf() const
	{
		return (length() == 1 && begin()->key.is_unity());
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline double BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseNorm(const ArgsTuple &argsTuple) const
	{
		const const_iterator itf = end();
		double retval = 0;
		
        for (const_iterator it = begin(); it != itf; ++it)
        {
			retval += it->cf.norm(argsTuple) * it->key.norm(argsTuple);
		}

		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::EvalType
		BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseEval(const double &t, const ArgsTuple &argsTuple) const
	{
		const const_iterator itf = end();
		EvalType retval(0);

		for (const_iterator it = begin(); it != itf; ++it)
        {
			EvalType tmp(it->cf.eval(t, argsTuple));
			tmp *= it->key.eval(t, argsTuple);
			retval += tmp;
		}

		return retval;
	}


	/// Length of the series.
	/**
	 * @return number of terms in the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::size_type BaseSeries<__PIRANHA_BASE_SERIES_TP>::length() const
	{
		return container.size();
	}


	/// Is series empty?
	/**
	 * @return true if the series has zero terms, false otherwise.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool BaseSeries<__PIRANHA_BASE_SERIES_TP>::empty() const
	{
		return container.empty();
	}


	/// Number of atoms in the series.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::size_type BaseSeries<__PIRANHA_BASE_SERIES_TP>::atoms() const
	{
		size_type retval = 0;
		const const_iterator itf = end();
		for (const_iterator it = begin(); it != itf; ++it)
        {
			retval += it->cf.atoms() + it->key.atoms();
		}

		return retval;
	}


	// TODO: update docs on equality test methods.
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
	template <class T>
	inline bool BaseSeries<__PIRANHA_BASE_SERIES_TP>::genericSeriesComparison(const T &other) const
	{
		if (length() != other.length())
        {
			return false;
		}

		const const_iterator itf = end();
        const const_iterator itfOther = other.end();
		
        for (const_iterator it = begin(); it != itf; ++it)
        {
			const_iterator itOther(other.findTerm(*it));
			if (itOther == itfOther || !(itOther->cf == it->cf))
            {
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
	inline bool BaseSeries<__PIRANHA_BASE_SERIES_TP>::genericNumericalComparison(const Number &x) const
	{
		if (numerical_comparison_zero_check<Number>()(x))
        {
			return empty();
		}

		if (!isSingleCf())
        {
			return false;
		}

		return (begin()->cf == x);
	}


	/// Base numerical equality test.
	/**
	 * @param[in] x argument of comparison.
	 *
	 * @return genericNumericalComparison() on x.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T>
	inline bool BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseEqualTo(const T &x) const
	{
		return BaseSeriesEqualToSelector<T>::run(*derived_const_cast, x);
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
