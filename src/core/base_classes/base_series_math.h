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


#include "../config.h"
#include "../exceptions.h"
#include "../math.h"
#include "../mp.h"
#include "../settings.h"
#include "../type_traits.h"
#include "base_series_def.h"
#include "base_series_mp.h"
#include "base_series_tag.h"

#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_complex.hpp>

#include <cstddef>
#include <complex>
#include <string>
#include <memory>
#include <type_traits>

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	// Do not use this to merge with self, assertion will fail.
    //merge terms of two series. wa precisely does that mean. make different types of terms from different series homogenous
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Derived2, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::mergeTerms(const Derived2 &series2, const ArgsTuple &argsTuple)
	{
		typedef typename Derived2::const_iterator const_iterator2;

		PIRANHA_ASSERT((void *)(std::addressof(*derived_cast)) != (void *)(std::addressof(series2))); // check on not insert into itself
		
        const const_iterator2 it_f = series2.end();
		for (const_iterator2 it = series2.begin(); it != it_f; ++it) 
		{
			// No need to check, we are merging from another series.
			insert<false, Sign>(*it, argsTuple); //without canonical check
		}

		return *derived_cast;
	}


	template <int N>
	struct MultDivCoefficientsHelper
	{
        static_assert(N == 0, "N must be either 0 or 1."); // this checks on 0 , 1 is differentn specialisation. see below.

		template <class Cf, class T, class ArgsTuple>
		static void run(Cf &cf, const T &x, const ArgsTuple &argsTuple) 
		{
			cf.multBy(x, argsTuple);
		}
	};


	template <>
	struct MultDivCoefficientsHelper<1>
	{
		template <class Cf, class T, class ArgsTuple>
		static void run(Cf &cf, const T &x, const ArgsTuple &argsTuple)
        {
			cf.divideBy(x, argsTuple);
		}
	};

    //
	// Multiply (N = 0) or divide (N = 1) all the coefficients of the series by a generic quantity x.
    // handle multiplication adn division in one implementation
    //
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <int N, class T, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::multDivCoefficientsBy(const T &x, const ArgsTuple &argsTuple)
	{
		const const_iterator itf = end();
		bool needsRebuilding = false;
		const_iterator it = begin();
		for (; it != itf; ++it) 
		{
			MultDivCoefficientsHelper<N>::run(it->cf, x, argsTuple);
			// If a term becomes ignorable once its coefficient has been multiplied/divided,
			// set the needs_rebuilding flag to true, increase the iterator and break out.
			if (it->cf.isIgnorable(argsTuple)) 
			{
				needsRebuilding = true;
				++it;
				break;
			}
		}

		// In case we broke out the cycle above, take care of multiplying/dividing the remaining
		// terms without checking for ignorability.
		for (; it != itf; ++it) 
		{
			MultDivCoefficientsHelper<N>::run(it->cf, x, argsTuple);
		}

		// Rebuild if needed.
		if (needsRebuilding) 
		{
			// TODO: probably insert range here.
			Derived newSeries;
			for (it = begin(); it != itf; ++it) 
			{
				newSeries.insert(*it, argsTuple);
			}

			baseSwap(newSeries);
		}
	}


	/// Merge series with a number.
	/**
	 * Term is constructed from coefficient constructed from number and default key,
	 * and then inserted into the series.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class Number, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::mergeWithNumber(const Number &n, const ArgsTuple &argsTuple)
	{
		typename Derived::TermType term(typename Derived::TermType::CfType(n, argsTuple), typename Derived::TermType::KeyType());

		insert<true, Sign>(term, argsTuple);
		return *derived_cast;
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::baseAdd(const T &x, const ArgsTuple &argsTuple)
	{
		return BaseSeriesAddSelector<T>::run(*derived_cast, x, argsTuple);
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSubtract(const T &x, const ArgsTuple &argsTuple)
	{
		return BaseSeriesSubtractSelector<T>::run(*derived_cast, x, argsTuple);
	}


	template<typename T>
	concept ComplexNumerical =  boost::is_complex<T>::value && 
					(std::is_integral_v<typename T::value_type> ||
					std::is_floating_point_v<typename T::value_type>);

	// This helper is needed because the standard complex<> class is missing comparison
	// operators against scalar types different from complex<>::value_type - whereas this
	// comparison is instead available in both the mp classes and the numerical coefficient classes.
	template <typename T>
	struct MultDivCoefficientsChecker
	{
		static bool checkZero(const T &x)
		{
			return x == 0;
		}

		static bool checkNonUnitary(const T &x)
		{
			return x != 1;
		}
	};


	// Enable for complex C++ integral and floating point types.
	// NOTE: it might be possible to use generically the code below exclusively,
	// but in case T is a series, extraction of real and imaginary part would be
	// quite a bit more expensive. So this also serves as an optimisation.
	template <ComplexNumerical T>
	struct MultDivCoefficientsChecker<T>
	{
		static bool checkZero(const T &c)
		{
			return (c.real() == 0 && c.imag() == 0);
		}

		static bool checkNonUnitary(const T &c)
		{
			return (c.real() != 1 || c.imag() != 0);
		}
	};


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::multiplyCoefficientsBy(const Number &x,
		const ArgsTuple &argsTuple)
	{
		if (MultDivCoefficientsChecker<Number>::checkZero(x)) 
		{
			clearTerms();

		} else if (MultDivCoefficientsChecker<Number>::checkNonUnitary(x)) 
		{
			multDivCoefficientsBy<0>(x,argsTuple);
		}

		return *derived_cast;
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::baseMultBy(const T &x, const ArgsTuple &argsTuple)
	{
		return BaseSeriesMultiplySelector<Derived, T>::run(*derived_cast, x, argsTuple);
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::divideCoefficientsBy(const Number &x, const ArgsTuple &argsTuple)
	{
		if (MultDivCoefficientsChecker<Number>::checkZero(x)) 
		{
			PIRANHA_THROW(zero_division_error, "cannot divide by zero");

		} else if (MultDivCoefficientsChecker<Number>::checkNonUnitary(x)) 
		{
			multDivCoefficientsBy<1>(x, argsTuple); // 1 == division
		}

		return *derived_cast;
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &BaseSeries<PIRANHA_BASE_SERIES_TP>::baseDivideBy(const T &x, const ArgsTuple &argsTuple)
	{
        static_assert((!PiranhaSeries<T>), "Cannot divide by another series.");

		return divideCoefficientsBy(x, argsTuple);
	}


    // partial derivative
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Series, class PosTuple, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::basePartial(const Derived &in, Series &out, const PosTuple &posTuple, const ArgsTuple &argsTuple)
	{
        static_assert(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
			                 "Size mismatch between args tuple and pos tuple in partial derivative.");

		PIRANHA_ASSERT(out.empty());

		const const_iterator itf = in.end();

		for (const_iterator it = in.begin(); it != itf; ++it) 
		{
			Series tmp1(it->cf.template partial<Series>(posTuple, argsTuple));

			tmp1.baseMultBy(Series::baseSeriesFromKey(it->key, argsTuple), argsTuple);

			Series tmp2(it->key.template partial<Series>(posTuple, argsTuple));

			tmp2.baseMultBy(Series::baseSeriesFromCf(it->cf, argsTuple), argsTuple);

			out.baseAdd(tmp1, argsTuple);
			out.baseAdd(tmp2, argsTuple);
		}
	}


	/// n-th Partial derivative.
	/**
	 * Calls partial() on all terms of the series, and inserts the resulting terms into return value.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::basePartial(int n, const PosTuple &posTuple, const ArgsTuple &argsTuple) const
	{
		if (n < 0) 
		{
			PIRANHA_THROW(value_error, "for an n-th partial derivative, n must be non-negative");
		}

		Derived retval(*derived_const_cast);
		for (; n > 0; --n) 
		{
			Derived tmp;
			basePartial(retval, tmp, posTuple, argsTuple);
			tmp.baseSwap(retval);
		}

		return retval;
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::basePartial(const PosTuple &posTuple, const ArgsTuple &argsTuple) const
	{
		return basePartial(1, posTuple, argsTuple);
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline bool BaseSeries<PIRANHA_BASE_SERIES_TP>::commonPowHandler(const Number &y, Derived &retval, const ArgsTuple &argsTuple) const
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
			retval.insert(TermType(begin()->cf.pow(y, argsTuple), typename TermType::KeyType()), argsTuple);
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
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::basePow(const double &y,
		const ArgsTuple &argsTuple) const
	{
		Derived retval;
		if (!commonPowHandler(y, retval, argsTuple)) 
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


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::basePow(const mp_rational &q, const ArgsTuple &argsTuple) const
	{
		Derived retval;
		if (!commonPowHandler(q, retval, argsTuple))
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

						} catch (const not_implemented_error &)
                        {
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


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::realPower(const double, const ArgsTuple &) const
	{
		PIRANHA_THROW(not_implemented_error, "real power for this series type has not been implemented");
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::negativeIntegerPower(const int n, const ArgsTuple &) const
	{
		(void)n;
		PIRANHA_ASSERT(n < 0);

		PIRANHA_THROW(not_implemented_error, "negative integer power for this series type has not been implemented");
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::rationalPower(const mp_rational &, const ArgsTuple &) const
	{
		PIRANHA_THROW(not_implemented_error, "rational power for this series type has not been implemented");
	}


	// Exponentiation to natural number.
	//
	// Exponentiation by squaring is used internally.
	//
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::naturalPower(const std::size_t n, const ArgsTuple &argsTuple) const
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
			    while (i)
                {
				    if (i & 1) 
                    {
					    retval.baseMultBy(tmp, argsTuple);
					    --i;
				    }

				    i /= 2;
				    if (i != 0) 
                    {
					    tmp.baseMultBy(tmp, argsTuple);
				    }
			    }
			}
		}
	}

	return retval;
    }


	/// Nth root.
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::baseRoot(const int n, const ArgsTuple &argsTuple) const
	{
		return basePow(mp_rational(1, n), argsTuple);
	}


	// Series inversion will use exponentiation to -1 as default.
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::baseInvert(const ArgsTuple &) const
	{
		PIRANHA_THROW(not_implemented_error, "inversion for this series type has not been implemented");
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
