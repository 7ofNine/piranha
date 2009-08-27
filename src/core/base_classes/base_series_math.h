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
#include "../math.h"
#include "../mp.h"
#include "../settings.h"
#include "base_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	// Do not use this to merge with self, assertion will fail.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Derived2, class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::merge_terms(const Derived2 &s2, const ArgsTuple &args_tuple)
	{
		typedef typename Derived2::const_iterator const_iterator2;
		piranha_assert((void *)derived_cast != (void *)&s2);
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
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::mult_div_coefficients_by(const T &x,
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
			// TODO: probably insert range here.
			Derived new_series;
			for (it = begin(); it != it_f; ++it) {
				new_series.insert(*it,args_tuple);
			}
			base_swap(new_series);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::multiply_coefficients_by(const T &x,
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
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::merge_with_number(const Number &n,
		const ArgsTuple &args_tuple)
	{
		typename Derived::term_type term(typename Derived::term_type::cf_type(n, args_tuple),
			typename Derived::term_type::key_type());
		insert<true, Sign>(term, args_tuple);
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_add(const Derived &s2, const ArgsTuple &args_tuple)
	{
		return merge_terms<true>(s2, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_add(const double &x, const ArgsTuple &args_tuple)
	{
		return merge_with_number<true>(x, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_add(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return merge_with_number<true>(q, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_add(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return merge_with_number<true>(z, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_subtract(const double &x, const ArgsTuple &args_tuple)
	{
		return merge_with_number<false>(x, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_subtract(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return merge_with_number<false>(q, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_subtract(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return merge_with_number<false>(z, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_subtract(const Derived &s2, const ArgsTuple &args_tuple)
	{
		return merge_terms<false>(s2, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::multiply_by_number(const Number &x,
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
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_mult_by(const double &x, const ArgsTuple &args_tuple)
	{
		return multiply_by_number(x, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_mult_by(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return multiply_by_number(q, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_mult_by(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return multiply_by_number(z, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_mult_by(const Derived &s2, const ArgsTuple &args_tuple)
	{
		derived_cast->multiply_by_series(s2, args_tuple);
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::divide_coefficients_by(const T &x,
			const ArgsTuple &args_tuple) {
		mult_div_coefficients_by<1>(x,args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::divide_by_number(const Number &x,
		const ArgsTuple &args_tuple)
	{
		if (x == 1) {
			return *derived_cast;
		} else if (x == 0) {
			piranha_throw(zero_division_error,"cannot divide by zero");
		}
		divide_coefficients_by(x, args_tuple);
		return *derived_cast;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_divide_by(const double &x, const ArgsTuple &args_tuple)
	{
		return divide_by_number(x, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_divide_by(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return divide_by_number(q, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_divide_by(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return divide_by_number(z, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Series, class PosTuple, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::ll_partial(const Derived &in, Series &out,
		const PosTuple &pos_tuple, const ArgsTuple &args_tuple)
	{
		p_static_check(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
			"Size mismatch between args tuple and pos tuple in partial derivative.");
		piranha_assert(out.empty());
		term_type tmp_term1, tmp_term2;
		const const_iterator it_f = in.end();
		for (const_iterator it = in.begin(); it != it_f; ++it) {
			Series tmp1(it->m_cf.template partial<Series>(pos_tuple,args_tuple));
			tmp1.base_mult_by(Series::base_series_from_key(it->m_key,args_tuple),args_tuple);
			Series tmp2(it->m_key.template partial<Series>(pos_tuple,args_tuple));
			tmp2.base_mult_by(Series::base_series_from_cf(it->m_cf,args_tuple),args_tuple);
			out.base_add(tmp1,args_tuple);
			out.base_add(tmp2,args_tuple);
		}
	}

	/// Partial derivative.
	/**
	 * Calls partial() on all terms of the series, and inserts the resulting terms into return value.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_partial(int n,
		const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const
	{
		if (n < 0) {
			piranha_throw(value_error,"for an n-th partial derivative, n must be non-negative");
		}
		Derived retval(*derived_const_cast);
		for (; n > 0; --n) {
			Derived tmp;
			ll_partial(retval,tmp,pos_tuple,args_tuple);
			tmp.base_swap(retval);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_partial(const PosTuple &pos_tuple,
		const ArgsTuple &args_tuple) const
	{
		return base_partial(1,pos_tuple,args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline bool toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::common_pow_handler(const Number &y, Derived &retval, const ArgsTuple &args_tuple) const
	{
		piranha_assert(retval.empty());
		// Handle the case of an empty series.
		if (empty()) {
			if (y < 0) {
				piranha_throw(zero_division_error,"cannot raise to negative power an empty series");
			} else if (y == 0) {
				// 0**0 == 1.
				retval.base_add(1, args_tuple);
				return true;
			} else {
				// 0**n == 0, with n > 0.
				return true;
			}
		}
		// If the series is a single cf, let's try to forward the pow call to the only coefficient.
		if (is_single_cf()) {
			retval.insert(term_type(begin()->m_cf.pow(y,args_tuple),typename term_type::key_type()),args_tuple);
			return true;
		}
		if (length() == 1) {
			// If length is 1, let's try to compute separately pow of coefficient and key
			// and assemble them. If we are not able, some exception is raised, we
			// catch it and try to go to the next step.
			try {
				const const_iterator it = begin();
				retval.insert(term_type(it->m_cf.pow(y, args_tuple), it->m_key.pow(y, args_tuple)), args_tuple);
				return true;
			} catch (...) {}
		}
		return false;
	}

	/// Real exponentiation.
	// Internally it will check if y is a real or an integer number.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_pow(const double &y,
		const ArgsTuple &args_tuple) const
	{
		Derived retval;
		if (!common_pow_handler(y, retval, args_tuple)) {
			if (is_integer(y)) {
				const int n = (int)(y);
				// NOTE: check if it is possible to skip using temporary variables.
				// Apparently this will be possible in C++1x with the move constructor.
				if (n >= 0) {
					Derived tmp(derived_const_cast->natural_power((size_t)n, args_tuple));
					retval.base_swap(tmp);
				} else {
					if (n == -1) {
						// Take advantage of a re-implementation of inversion in derived class,
						// if available. Otherwise just use negative_integer_power.
						try {
							Derived tmp(derived_const_cast->base_inv(args_tuple));
							retval.base_swap(tmp);
						} catch (const not_implemented_error &) {
							Derived tmp(derived_const_cast->negative_integer_power(n, args_tuple));
							retval.base_swap(tmp);
						}
					} else {
						Derived tmp(derived_const_cast->negative_integer_power(n, args_tuple));
						retval.base_swap(tmp);
					}
				}
			} else {
				Derived tmp(derived_const_cast->real_power(y, args_tuple));
				retval.base_swap(tmp);
			}
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_pow(const mp_rational &q,
		const ArgsTuple &args_tuple) const
	{
		Derived retval;
		if (!common_pow_handler(q, retval, args_tuple)) {
			// If rational is integer, dispatch to natural power or negative integer power.
			if (q.get_den() == 1) {
				const int n = q.get_num().to_int();
				if (n >= 0) {
					Derived tmp(derived_const_cast->natural_power((size_t)n, args_tuple));
					retval.base_swap(tmp);
				} else {
					if (n == -1) {
						try {
							Derived tmp(derived_const_cast->base_inv(args_tuple));
							retval.base_swap(tmp);
						} catch (const not_implemented_error &) {
							Derived tmp(derived_const_cast->negative_integer_power(n, args_tuple));
							retval.base_swap(tmp);
						}
					} else {
						Derived tmp(derived_const_cast->negative_integer_power(n, args_tuple));
						retval.base_swap(tmp);
					}
				}
			} else {
				Derived tmp(derived_const_cast->rational_power(q, args_tuple));
				retval.base_swap(tmp);
			}
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::real_power(const double &, const ArgsTuple &) const
	{
		piranha_throw(not_implemented_error,"real power for this series type has not been implemented");
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::negative_integer_power(const int &n,
		const ArgsTuple &) const
	{
		(void)n;
		piranha_assert(n < 0);
		piranha_throw(not_implemented_error,"negative integer power for this series type has not been implemented");
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::rational_power(const mp_rational &, const ArgsTuple &) const
	{
		piranha_throw(not_implemented_error,"rational power for this series type has not been implemented");
	}

	/// Exponentiation to natural number.
	/**
	 * Exponentiation by squaring is used internally.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::natural_power(const size_t &n,
		const ArgsTuple &args_tuple) const
	{
		Derived retval;
		switch (n) {
		case 0: {
			retval.base_add(1, args_tuple);
			break;
		}
		case 1: {
			retval = *derived_const_cast;
			break;
		}
		case 2: {
			retval = *derived_const_cast;
			retval.base_mult_by(*derived_const_cast, args_tuple);
			break;
		}
		case 3: {
			retval = *derived_const_cast;
			retval.base_mult_by(*derived_const_cast, args_tuple);
			retval.base_mult_by(*derived_const_cast, args_tuple);
			break;
		}
		case 4: {
			retval = *derived_const_cast;
			retval.base_mult_by(*derived_const_cast, args_tuple);
			retval.base_mult_by(*derived_const_cast, args_tuple);
			retval.base_mult_by(*derived_const_cast, args_tuple);
			break;
		}
		default: {
			retval.base_add(1, args_tuple);
			// Use scoping here to have tmp destroyed when it is not needed anymore.
			{
			Derived tmp(*derived_const_cast);
			size_t i = n;
			while (i) {
				if (i & 1) {
					retval.base_mult_by(tmp, args_tuple);
					--i;
				}
				i /= 2;
				if (i != 0) {
					tmp.base_mult_by(tmp, args_tuple);
				}
			}
			}
		}
		}
		return retval;
	}

	/// Nth root.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_root(const int &n,
			const ArgsTuple &args_tuple) const
	{
		return base_pow(mp_rational(1,n),args_tuple);
	}

	// Series inversion will use exponentiation to -1 as default.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_inv(const ArgsTuple &) const
	{
		piranha_throw(not_implemented_error,"inversion for this series type has not been implemented");
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
