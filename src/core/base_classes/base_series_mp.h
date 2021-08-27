/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redis\bute it and/or modify  *
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

#ifndef PIRANHA_BASE_SERIES_MP_H
#define PIRANHA_BASE_SERIES_MP_H

#include "../config.h"
#include "../exceptions.h"
#include "../type_traits.h"
#include "base_series_tag.h"

#include <complex>
#include <vector>
#include <type_traits>


namespace piranha
{
	// Implementation of factory methods from coefficients and keys.
	template <class RequestedKey, class SeriesKey>
	struct SeriesFromKeyImpl
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedKey &key, const ArgsTuple &argsTuple)
		{
			typedef typename Series::TermType::CfType CfSeriesType;
			CfSeriesType cfSeries;
			SeriesFromKeyImpl<RequestedKey, typename CfSeriesType::TermType::KeyType>::run(cfSeries, key, argsTuple);

			retval.insert(typename Series::TermType(cfSeries, typename Series::TermType::KeyType()), argsTuple);
		}
	};


	template <class Key>
	struct SeriesFromKeyImpl<Key, Key>
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Key &key, const ArgsTuple &argsTuple)
		{
			retval.insert(typename Series::TermType(typename Series::TermType::CfType(1, argsTuple), key), argsTuple);
		}
	};


	template <class RequestedCf, class SeriesCf>
	struct SeriesFromCfImpl
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedCf &cf, const ArgsTuple &argsTuple)
		{
			typedef typename Series::TermType::CfType CfSeriesType;
			CfSeriesType cfSeries;
			SeriesFromCfImpl<RequestedCf, typename CfSeriesType::TermType::CfType>::run(cfSeries, cf, argsTuple);

			retval.insert(typename Series::TermType(cfSeries, typename Series::TermType::KeyType()), argsTuple);
		}
	};


	template <class Cf>
	struct SeriesFromCfImpl<Cf, Cf>
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Cf &cf, const ArgsTuple &argsTuple)
		{
			retval.insert(typename Series::TermType(cf, typename Series::TermType::KeyType()), argsTuple);
		}
	};


	// This functor will disentangle and build the flattened terms iterating recursively through the echelon levels.
	// N is the echelon level we are working on
	template <int N>
	class SeriesFlattener
	{
		public:
            static_assert(N > 0, "");

			template <class CfSeries, class Term, class ArgsTuple>
			static void run(CfSeries &cfSeries, Term &term, std::vector<Term> &out, const ArgsTuple &argsTuple)
			{
				// For each coefficient series (which is residing inside the insertion term), we create a copy of it,
				// then we insert one by one its original terms and, step by step, we go deeper into the recursion.
				PIRANHA_ASSERT(!cfSeries.empty());

				const CfSeries tmp(cfSeries);

				for (typename CfSeries::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
                {
					cfSeries.clearTerms();
					cfSeries.insert(*it, argsTuple);
					SeriesFlattener<N - 1>::run(cfSeries.begin()->cf, term, out, argsTuple);
				}
			}
	};


	template <>
	class SeriesFlattener<0>
	{
		public:
			template <class Cf, class Term, class ArgsTuple>
			static void run(Cf &, Term &term, std::vector<Term> &out, const ArgsTuple &)
			{
				out.push_back(term);
			}
	};


	// These structs are used to select at compile time which low-level methods in BaseSeries
	// to call to implement arithmetic operations - based on the type of argument.


	template<typename T>
	struct BaseSeriesAddSelector
	{
		template <typename Derived, typename ArgsTuple>
		static Derived &run(Derived &series, const T &x, const ArgsTuple &argsTuple)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.template mergeTerms<true>(x, argsTuple); // with sign
			}
			else
			{
				return series.template mergeWithNumber<true>(x, argsTuple);
			}
		}
	};


	template <typename T>
	struct BaseSeriesSubtractSelector
	{
		template <class Derived, class ArgsTuple>
		static Derived &run(Derived& series, const T& x, const ArgsTuple& argsTuple)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.template mergeTerms<false>(x, argsTuple); // subtract
			}
			else
			{
				return series.template mergeWithNumber<false>(x, argsTuple); // subtract 
			}
		}
	};


	template <typename Derived, typename T>
	struct BaseSeriesMultiplySelector
	{
		template <class ArgsTuple>
		static Derived &run(Derived &series, const T &x, const ArgsTuple &argsTuple)
		{
			return series.multiplyCoefficientsBy(x, argsTuple);
		}
	};


	template <typename T, typename  Derived>  requires std::is_base_of_v<BaseSeriesTag, T> && (std::is_same_v<Derived, T> || std::is_same_v<Derived, std::complex<T> >)
	struct BaseSeriesMultiplySelector<Derived, T>
	{
		template <typename ArgsTuple>
		static Derived &run(Derived &series, const T &other, const ArgsTuple &argsTuple)
		{
			//std::cout << "BaseSeriesMultiplication_selector::run : 0" << std::endl << std::flush;
			series.multiply_by_series(other, argsTuple);
			//std::cout << "BaseSeriesMultiplication_selector::run : 1" << std::endl << std::flush;
			return series;
		}
	};


    //////////////////////////////////////////////////////////////
    //
    //

	template <typename T>
	struct BaseSeriesEqualToSelector 
	{
		template <class Derived>
		static bool run(const Derived& series, const T& x)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.genericSeriesComparison(x);
			}
			else
			{
				return series.genericNumericalComparison(x);
			}
		}
	};

}

#endif
