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

#include <boost/lexical_cast.hpp>
#include <string>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../p_assert.h"
#include "../settings.h"

namespace piranha
{
	// Do not use this to merge with self, assertion will fail.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Derived2, class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::merge_terms(const Derived2 &s2, const ArgsTuple &args_tuple)
	{
		typedef typename Derived2::const_iterator const_iterator2;
		p_assert((void *)derived_cast != (void *)&s2);
		const const_iterator2 it_f = s2.end();
		for (const_iterator2 it = s2.begin(); it != it_f; ++it) {
			// No need to check, we are merging from another series.
			insert<false, Sign>(*it, args_tuple);
		}
		return *derived_cast;
	}

	template <int N>
	struct mult_div_coefficients_helper
	{
		p_static_check(N == 0, "N must be either 0 or 1.");
		template <class Cf, class T, class ArgsTuple>
		static void run(Cf &cf, const T &x, const ArgsTuple &args_tuple) {
			cf.mult_by(x,args_tuple);
		}
	};

	template <>
	struct mult_div_coefficients_helper<1>
	{
		template <class Cf, class T, class ArgsTuple>
		static void run(Cf &cf, const T &x, const ArgsTuple &args_tuple) {
			cf.divide_by(x,args_tuple);
		}
	};

	// Multiply (N = 0) or divide (N = 1) all the coefficients of the series by a generic quantity x.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <int N, class T, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::mult_div_coefficients_by(const T &x,
			const ArgsTuple &args_tuple)
	{
		const const_iterator it_f = end();
		bool needs_rebuilding = false;
		const_iterator it = begin();
		for (; it != it_f; ++it) {
			mult_div_coefficients_helper<N>::run(it->m_cf,x,args_tuple);
			// If a term becomes ignorable once its coefficient has been multiplied/divided,
			// set the needs_rebuilding flag to true, increase the iterator and break out.
			if (it->m_cf.is_ignorable(args_tuple)) {
				needs_rebuilding = true;
				++it;
				break;
			}
		}
		// In case we broke out the cycle above, take care of multiplying/dividing the remaining
		// terms without checking for ignorability.
		for (; it != it_f; ++it) {
			mult_div_coefficients_helper<N>::run(it->m_cf,x,args_tuple);
		}
		// Rebuild if needed.
		if (needs_rebuilding) {
			Derived new_series;
			for (it = begin(); it != it_f; ++it) {
				new_series.insert(*it,args_tuple);
			}
			swap_terms(new_series);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::multiply_coefficients_by(const T &x,
			const ArgsTuple &args_tuple) {
		mult_div_coefficients_by<0>(x,args_tuple);
	}

	/// Merge series with a number.
	/**
	 * Term is constructed from coefficient constructed from number and default key,
	 * and then inserted into the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Number, class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::merge_with_number(const Number &n,
		const ArgsTuple &args_tuple)
	{
		typename Derived::term_type term(typename Derived::term_type::cf_type(n, args_tuple),
			typename Derived::term_type::key_type());
		insert<true, Sign>(term, args_tuple);
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::add(const Derived &s2, const ArgsTuple &args_tuple)
	{
		return merge_terms<true>(s2, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::subtract(const Derived &s2, const ArgsTuple &args_tuple)
	{
		return merge_terms<false>(s2, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::multiply_by_number(const Number &x,
		const ArgsTuple &args_tuple)
	{
		if (x == 0) {
			clear_terms();
		} else if (x != 1) {
			multiply_coefficients_by(x, args_tuple);
		}
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::mult_by(const max_fast_int &n, const ArgsTuple &args_tuple)
	{
		return multiply_by_number(n, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::mult_by(const double &x, const ArgsTuple &args_tuple)
	{
		return multiply_by_number(x, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::mult_by(const Derived &s2, const ArgsTuple &args_tuple)
	{
		derived_cast->multiply_by_series(s2, args_tuple);
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::divide_coefficients_by(const T &x,
			const ArgsTuple &args_tuple) {
		mult_div_coefficients_by<1>(x,args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::divide_by_number(const Number &x,
		const ArgsTuple &args_tuple)
	{
		if (x == 1) {
			return *derived_cast;
		} else if (x == 0) {
			throw division_by_zero();
		}
		divide_coefficients_by(x, args_tuple);
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::divide_by(const max_fast_int &n, const ArgsTuple &args_tuple)
	{
		return divide_by_number(n, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &base_series<__PIRANHA_BASE_SERIES_TP>::divide_by(const double &x, const ArgsTuple &args_tuple)
	{
		return divide_by_number(x, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::choose(
		const max_fast_int &n_, const max_fast_int &k, const ArgsTuple &args_tuple)
	{
		const max_fast_int n = (n_ >= 0) ? n_ : k - n_ - 1;
		if (k < 0 || k > n) {
			throw unsuitable("Invalid k for binomial coefficient.");
		}
		// Starting point is 1.
		Derived retval;
		retval.insert(term_type(cf_type(max_fast_int(1),args_tuple),key_type()),args_tuple);
		for (max_fast_int i = 1; i <= k; ++i) {
			retval.mult_by(n - k + i, args_tuple);
			retval.divide_by(i, args_tuple);
		}
		if (n_ < 0) {
			retval.mult_by(cs_phase(k), args_tuple);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::factorial(const max_fast_int &n, const ArgsTuple &args_tuple)
	{
		if (n < 0) {
			throw unsuitable("Negative argument for factorial.");
		}
		Derived retval;
		retval.insert(term_type(cf_type(max_fast_int(1),args_tuple),key_type()),args_tuple);
		for (max_fast_int i = 2; i <= n; ++i) {
			retval.mult_by(i,args_tuple);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::ll_partial(const Derived &in, Derived &out,
		const PosTuple &pos_tuple, const ArgsTuple &args_tuple) {
		p_static_check(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
			"Size mismatch between args tuple and pos tuple in partial derivative.");
		p_assert(out.empty());
		term_type tmp_term1, tmp_term2;
		const const_iterator it_f = in.end();
		for (const_iterator it = in.begin(); it != it_f; ++it) {
			it->partial(tmp_term1, tmp_term2, pos_tuple, args_tuple);
			out.insert(tmp_term1, args_tuple);
			out.insert(tmp_term2, args_tuple);
		}
	}

	/// Partial derivative.
	/**
	 * Calls partial() on all terms of the series, and inserts the resulting terms into return value.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::partial(max_fast_int n,
		const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const
	{
		if (n < 0) {
			throw unsuitable("For an n-th partial derivative, n must be non-negative.");
		}
		Derived retval(*derived_const_cast);
		for (; n > 0; --n) {
			Derived tmp;
			ll_partial(retval,tmp,pos_tuple,args_tuple);
			tmp.swap_terms(retval);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::partial(const PosTuple &pos_tuple,
		const ArgsTuple &args_tuple) const {
		return partial(1,pos_tuple,args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline bool base_series<__PIRANHA_BASE_SERIES_TP>::common_power_handler(const Number &y, Derived &retval,
			const ArgsTuple &args_tuple) const
	{
		p_assert(retval.empty());
		// Handle the case of an empty series.
		if (empty()) {
			if (y < 0) {
				throw division_by_zero();
			} else if (y == 0) {
				// 0**0 == 1.
				retval = Derived((max_fast_int)1, args_tuple);
				return true;
			} else {
				// 0**n == 0, with n > 0.
				return true;
			}
		} else if (length() == 1) {
			// If the series has a single term, dispatch pow to the coefficient and key of said term.
			const const_iterator it = begin();
			retval.insert(term_type(it->m_cf.pow(y, args_tuple), it->m_key.pow(y, args_tuple)), args_tuple);
			return true;
		}
		return false;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::pow(const max_fast_int &n,
			const ArgsTuple &args_tuple) const
	{
		Derived retval;
		if (n >= 0) {
			Derived tmp(derived_const_cast->natural_power((size_t)n, args_tuple));
			retval.swap_terms(tmp);
		} else if (!common_power_handler(n, retval, args_tuple)) {
			Derived tmp(derived_const_cast->negative_integer_power(n, args_tuple));
			retval.swap_terms(tmp);
		}
		return retval;
	}

	/// Real exponentiation.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::pow(const double &y,
			const ArgsTuple &args_tuple) const
	{
		Derived retval;
		if (!common_power_handler(y, retval, args_tuple)) {
			Derived tmp(derived_const_cast->real_power(y, args_tuple));
			retval.swap_terms(tmp);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::real_power(const double &, const ArgsTuple &) const
	{
		throw(not_implemented("Real power for this series has not been implemented."));
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::negative_integer_power(const max_fast_int &n,
		const ArgsTuple &) const
	{
		(void)n;
		p_assert(n < 0);
		throw(not_implemented("Negative integer power for this series has not been implemented."));
	}

	/// Exponentiation to natural number.
	/**
	 * Exponentiation by squaring is used internally.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::natural_power(const size_t &n, const ArgsTuple &args_tuple) const
	{
		Derived retval;
		switch (n) {
		case 0: {
			retval = Derived(static_cast<max_fast_int>(1), args_tuple);
			break;
		}
		case 1: {
			retval = *derived_const_cast;
			break;
		}
		case 2: {
			retval = *derived_const_cast;
			retval.mult_by(*derived_const_cast, args_tuple);
			break;
		}
		case 3: {
			retval = *derived_const_cast;
			retval.mult_by(*derived_const_cast, args_tuple);
			retval.mult_by(*derived_const_cast, args_tuple);
			break;
		}
		case 4: {
			retval = *derived_const_cast;
			retval.mult_by(*derived_const_cast, args_tuple);
			retval.mult_by(*derived_const_cast, args_tuple);
			retval.mult_by(*derived_const_cast, args_tuple);
			break;
		}
		default: {
			retval = Derived(static_cast<max_fast_int>(1), args_tuple);
			// Use scoping here to have tmp destroyed when it is not needed anymore.
			{
			Derived tmp(*derived_const_cast);
			size_t i = n;
			while (i) {
				if (i & 1) {
					retval.mult_by(tmp, args_tuple);
					--i;
				}
				i /= 2;
				if (i != 0) {
					tmp.mult_by(tmp, args_tuple);
				}
			}
			}
		}
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool base_series<__PIRANHA_BASE_SERIES_TP>::common_root_handler(const max_fast_int &n, Derived &retval,
			const ArgsTuple &args_tuple) const
	{
		p_assert(retval.empty());
		if (n == 0) {
			throw division_by_zero();
		}
		if (n == 1) {
			retval = *derived_const_cast;
			return true;
		}
		// Handle the case of an empty series.
		if (empty()) {
			if (n < 0) {
				throw division_by_zero();
			} else {
				// 0**n == 0, with n > 0.
				return true;
			}
		} else if (length() == 1) {
			// If the series has a single term, dispatch pow to the coefficient and key of said term.
			const const_iterator it = begin();
			retval.insert(term_type(it->m_cf.root(n, args_tuple), it->m_key.root(n, args_tuple)), args_tuple);
			return true;
		}
		return false;
	}

	/// Nth root.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::root(const max_fast_int &n,
			const ArgsTuple &args_tuple) const
	{
		Derived retval;
		if (!common_root_handler(n, retval, args_tuple)) {
			Derived tmp(derived_const_cast->nth_root(n, args_tuple));
			retval.swap_terms(tmp);
		}
		return retval;
	}

	// By default use real power.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::nth_root(const max_fast_int &n,
			const ArgsTuple &args_tuple) const
	{
		p_assert(n > 0);
		return pow(1. / static_cast<double>(n), args_tuple);
	}
}

#endif
