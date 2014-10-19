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
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/addressof.hpp>
#include <boost/utility/enable_if.hpp>
#include <cstddef>
#include <complex>
#include <string>

#include "../config.h"
#include "../exceptions.h"
#include "../math.h"
#include "../mp.h"
#include "../settings.h"
#include "base_series_def.h"
#include "base_series_mp.h"
#include "base_series_tag.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	// Do not use this to merge with self, assertion will fail.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Derived2, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::mergeTerms(const Derived2 &s2, const ArgsTuple &argsTuple)
	{
		typedef typename Derived2::const_iterator const_iterator2;
		PIRANHA_ASSERT((void *)boost::addressof(*derived_cast) != (void *)boost::addressof(s2));
		const const_iterator2 it_f = s2.end();
		for (const_iterator2 it = s2.begin(); it != it_f; ++it) 
		{
			// No need to check, we are merging from another series.
			insert<false, Sign>(*it, argsTuple);
		}
		return *derived_cast;
	}


	template <int N>
	struct mult_div_coefficients_helper
	{
		PIRANHA_STATIC_CHECK(N == 0, "N must be either 0 or 1.");
		template <class Cf, class T, class ArgsTuple>
		static void run(Cf &cf, const T &x, const ArgsTuple &argsTuple) 
		{
			cf.mult_by(x,argsTuple);
		}
	};


	template <>
	struct mult_div_coefficients_helper<1>
	{
		template <class Cf, class T, class ArgsTuple>
		static void run(Cf &cf, const T &x, const ArgsTuple &argsTuple)
        {
			cf.divideBy(x,argsTuple);
		}
	};


	// Multiply (N = 0) or divide (N = 1) all the coefficients of the series by a generic quantity x.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <int N, class T, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::multDivCoefficientsBy(const T &x,
			const ArgsTuple &argsTuple)
	{
		const const_iterator itf = end();
		bool needs_rebuilding = false;
		const_iterator it = begin();
		for (; it != itf; ++it) 
		{
			mult_div_coefficients_helper<N>::run(it->cf, x, argsTuple);
			// If a term becomes ignorable once its coefficient has been multiplied/divided,
			// set the needs_rebuilding flag to true, increase the iterator and break out.
			if (it->cf.is_ignorable(argsTuple)) 
			{
				needs_rebuilding = true;
				++it;
				break;
			}
		}

		// In case we broke out the cycle above, take care of multiplying/dividing the remaining
		// terms without checking for ignorability.
		for (; it != itf; ++it) 
		{
			mult_div_coefficients_helper<N>::run(it->cf, x, argsTuple);
		}

		// Rebuild if needed.
		if (needs_rebuilding) 
		{
			// TODO: probably insert range here.
			Derived new_series;
			for (it = begin(); it != itf; ++it) 
			{
				new_series.insert(*it, argsTuple);
			}

			baseSwap(new_series);
		}
	}


	/// Merge series with a number.
	/**
	 * Term is constructed from coefficient constructed from number and default key,
	 * and then inserted into the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Number, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::mergeWithNumber(const Number &n, const ArgsTuple &argsTuple)
	{
		typename Derived::TermType term(typename Derived::TermType::cf_type(n, argsTuple), typename Derived::TermType::key_type());

		insert<true, Sign>(term, argsTuple);
		return *derived_cast;
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseAdd(const T &x, const ArgsTuple &argsTuple)
	{
		return BaseSeriesAddSelector<T>::run(*derived_cast, x, argsTuple);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSubtract(const T &x, const ArgsTuple &argsTuple)
	{
		return BaseSeriesSubtractSelector<T>::run(*derived_cast, x, argsTuple);
	}

	// This helper is needed because the standard complex<> class is missing comparison
	// operators against scalar types different from complex<>::value_type - whereas this
	// comparison is instead available in both the mp classes and the numerical coefficient classes.
	template <class T, class Enable = void>
	struct multdiv_coefficients_helper
	{
		static bool check_zero(const T &x)
		{
			return x == 0;
		}
		static bool check_non_unitary(const T &x)
		{
			return x != 1;
		}
	};


	// Enable for complex C++ integral and floating point types.
	// NOTE: it might be possible to use generically the code below exclusively,
	// but in case T is a series, extraction of real and imaginary part would be
	// quite a bit more expensive. So this also serves as an optimisation.
	template <class T>
	struct multdiv_coefficients_helper<T,typename boost::enable_if_c<boost::is_complex<T>::value && (
		boost::is_integral<typename T::value_type>::value ||
		boost::is_floating_point<typename T::value_type>::value)>::type>
	{
		static bool check_zero(const T &c)
		{
			return (c.real() == 0 && c.imag() == 0);
		}

		static bool check_non_unitary(const T &c)
		{
			return (c.real() != 1 || c.imag() != 0);
		}
	};


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::multiplyCoefficientsBy(const Number &x,
		const ArgsTuple &argsTuple)
	{
		if (multdiv_coefficients_helper<Number>::check_zero(x)) 
		{
			clearTerms();
		} else if (multdiv_coefficients_helper<Number>::check_non_unitary(x)) 
		{
			multDivCoefficientsBy<0>(x,argsTuple);
		}
		return *derived_cast;
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseMultBy(const T &x, const ArgsTuple &argsTuple)
	{
		return BaseSeriesMultiplySelector<Derived, T>::run(*derived_cast, x, argsTuple);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::divideCoefficientsBy(const Number &x,
		const ArgsTuple &argsTuple)
	{
		if (multdiv_coefficients_helper<Number>::check_zero(x)) 
		{
			PIRANHA_THROW(zero_division_error,"cannot divide by zero");

		} else if (multdiv_coefficients_helper<Number>::check_non_unitary(x)) 
		{
			multDivCoefficientsBy<1>(x, argsTuple);
		}
		return *derived_cast;
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseDivideBy(const T &x, const ArgsTuple &argsTuple)
	{
		PIRANHA_STATIC_CHECK((!boost::is_base_of<BaseSeriesTag,T>::value),"Cannot divide by another series.");

		return divideCoefficientsBy(x, argsTuple);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Series, class PosTuple, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::basePartial(const Derived &in, Series &out, const PosTuple &pos_tuple, const ArgsTuple &argsTuple)
	{
		PIRANHA_STATIC_CHECK(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
			                 "Size mismatch between args tuple and pos tuple in partial derivative.");
		PIRANHA_ASSERT(out.empty());

		TermType tmp_term1;
        TermType tmp_term2;
		const const_iterator itf = in.end();

		for (const_iterator it = in.begin(); it != itf; ++it) 
		{
			Series tmp1(it->cf.template partial<Series>(pos_tuple, argsTuple));

			tmp1.baseMultBy(Series::baseSeriesFromKey(it->key, argsTuple), argsTuple);

			Series tmp2(it->key.template partial<Series>(pos_tuple, argsTuple));

			tmp2.baseMultBy(Series::baseSeriesFromCf(it->cf, argsTuple), argsTuple);

			out.baseAdd(tmp1, argsTuple);
			out.baseAdd(tmp2, argsTuple);
		}
	}


	/// Partial derivative.
	/**
	 * Calls partial() on all terms of the series, and inserts the resulting terms into return value.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::basePartial(int n, const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
	{
		if (n < 0) 
		{
			PIRANHA_THROW(value_error,"for an n-th partial derivative, n must be non-negative");
		}

		Derived retval(*derived_const_cast);
		for (; n > 0; --n) 
		{
			Derived tmp;
			basePartial(retval, tmp, pos_tuple, argsTuple);
			tmp.baseSwap(retval);
		}

		return retval;
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::basePartial(const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
	{
		return basePartial(1, pos_tuple, argsTuple);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline bool BaseSeries<__PIRANHA_BASE_SERIES_TP>::common_pow_handler(const Number &y, Derived &retval, const ArgsTuple &argsTuple) const
	{
		PIRANHA_ASSERT(retval.empty());
		// Handle the case of an empty series.
		if (empty()) 
		{
			if (y < 0) 
			{
				PIRANHA_THROW(zero_division_error,"cannot raise to negative power an empty series");

			} else if (y == 0) 
			{
				// 0**0 == 1.
				retval.baseAdd(1, argsTuple);
				return true;

			} else 
			{
				// 0**n == 0, with n > 0.
				return true;
			}
		}

		// If the series is a single cf, let's try to forward the pow call to the only coefficient.
		if (isSingleCf()) 
		{
			retval.insert(TermType(begin()->cf.pow(y, argsTuple), typename TermType::key_type()), argsTuple);
			return true;
		}
		if (length() == 1) 
		{
			// If length is 1, let's try to compute separately pow of coefficient and key
			// and assemble them. If we are not able, some exception is raised, we
			// catch it and try to go to the next step.
			try {
				const const_iterator it = begin();
				retval.insert(TermType(it->cf.pow(y, argsTuple), it->key.pow(y, argsTuple)), argsTuple);
				return true;

			} catch (...) {}
		}

		return false;
	}

	/// Real exponentiation.
	// Internally it will check if y is a real or an integer number.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::basePow(const double &y,
		const ArgsTuple &argsTuple) const
	{
		Derived retval;
		if (!common_pow_handler(y, retval, argsTuple)) 
		{
			if (is_integer(y)) 
			{
				const int n = (int)(y);
				// NOTE: check if it is possible to skip using temporary variables.
				// Apparently this will be possible in C++1x with the move constructor.
				if (n >= 0) 
				{
					Derived tmp(derived_const_cast->naturalPower((std::size_t)n, argsTuple));
					retval.baseSwap(tmp);

				} else 
				{
					if (n == -1) 
					{
						// Take advantage of a re-implementation of inversion in derived class,
						// if available. Otherwise just use negativeIntegerPower.
						try {
							Derived tmp(derived_const_cast->baseInvert(argsTuple));
							retval.baseSwap(tmp);

						} catch (const not_implemented_error &) {
							Derived tmp(derived_const_cast->negativeIntegerPower(n, argsTuple));
							retval.baseSwap(tmp);
						}
					} else 
					{
						Derived tmp(derived_const_cast->negativeIntegerPower(n, argsTuple));
						retval.baseSwap(tmp);
					}
				}
			} else 
			{
				Derived tmp(derived_const_cast->realPower(y, argsTuple));
				retval.baseSwap(tmp);
			}
		}

		return retval;
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::basePow(const mp_rational &q,
		const ArgsTuple &argsTuple) const
	{
		Derived retval;
		if (!common_pow_handler(q, retval, argsTuple))
		{
			// If rational is integer, dispatch to natural power or negative integer power.
			if (q.get_den() == 1) 
			{
				const int n = q.get_num().to_int();
				if (n >= 0) 
				{
					Derived tmp(derived_const_cast->naturalPower((std::size_t)n, argsTuple));
					retval.baseSwap(tmp);
				} else 
				{
					if (n == -1) 
					{
						try {
							Derived tmp(derived_const_cast->baseInvert(argsTuple));
							retval.baseSwap(tmp);
						} catch (const not_implemented_error &) {
							Derived tmp(derived_const_cast->negativeIntegerPower(n, argsTuple));
							retval.baseSwap(tmp);
						}
					} else 
					{
						Derived tmp(derived_const_cast->negativeIntegerPower(n, argsTuple));
						retval.baseSwap(tmp);
					}
				}
			} else 
			{
				Derived tmp(derived_const_cast->rationalPower(q, argsTuple));
				retval.baseSwap(tmp);
			}
		}
		return retval;
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::realPower(const double &, const ArgsTuple &) const
	{
		PIRANHA_THROW(not_implemented_error,"real power for this series type has not been implemented");
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::negativeIntegerPower(const int &n, const ArgsTuple &) const
	{
		(void)n;
		PIRANHA_ASSERT(n < 0);
		PIRANHA_THROW(not_implemented_error,"negative integer power for this series type has not been implemented");
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::rationalPower(const mp_rational &, const ArgsTuple &) const
	{
		PIRANHA_THROW(not_implemented_error,"rational power for this series type has not been implemented");
	}


	/// Exponentiation to natural number.
	/**
	 * Exponentiation by squaring is used internally.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::naturalPower(const std::size_t &n,
		const ArgsTuple &argsTuple) const
	{
		Derived retval;
		switch (n) {
		case 0: {
			retval.baseAdd(1, argsTuple);
			break;
		}
		case 1: {
			retval = *derived_const_cast;
			break;
		}
		case 2: {
			retval = *derived_const_cast;
			retval.baseMultBy(*derived_const_cast, argsTuple);
			break;
		}
		case 3: {
			retval = *derived_const_cast;
			retval.baseMultBy(*derived_const_cast, argsTuple);
			retval.baseMultBy(*derived_const_cast, argsTuple);
			break;
		}
		case 4: {
			retval = *derived_const_cast;
			retval.baseMultBy(*derived_const_cast, argsTuple);
			retval.baseMultBy(*derived_const_cast, argsTuple);
			retval.baseMultBy(*derived_const_cast, argsTuple);
			break;
		}
		default: {
			retval.baseAdd(1, argsTuple);
			// Use scoping here to have tmp destroyed when it is not needed anymore.
			{
			Derived tmp(*derived_const_cast);
			std::size_t i = n;
			while (i) {
				if (i & 1) {
					retval.baseMultBy(tmp, argsTuple);
					--i;
				}
				i /= 2;
				if (i != 0) {
					tmp.baseMultBy(tmp, argsTuple);
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
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseRoot(const int &n,
			const ArgsTuple &argsTuple) const
	{
		return basePow(mp_rational(1, n), argsTuple);
	}


	// Series inversion will use exponentiation to -1 as default.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseInvert(const ArgsTuple &) const
	{
		PIRANHA_THROW(not_implemented_error,"inversion for this series type has not been implemented");
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
