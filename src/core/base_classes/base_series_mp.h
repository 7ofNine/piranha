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

#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "base_series_tag.h"

namespace piranha
{
	// Implementation of factory methods from coefficients and keys.
	template <class RequestedKey, class SeriesKey>
	struct SeriesFromKeyImpl
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedKey &key, const ArgsTuple &argsTuple)
		{
			typedef typename Series::TermType::cf_type cf_series_type;
			cf_series_type cf_series;
			SeriesFromKeyImpl<RequestedKey,typename cf_series_type::TermType::key_type>::run(cf_series, key, argsTuple);

			retval.insert(typename Series::TermType(cf_series, typename Series::TermType::key_type()), argsTuple);
		}
	};


	template <class Key>
	struct SeriesFromKeyImpl<Key,Key>
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Key &key, const ArgsTuple &argsTuple)
		{
			retval.insert(typename Series::TermType(typename Series::TermType::cf_type(1, argsTuple), key), argsTuple);
		}
	};


	template <class RequestedCf, class SeriesCf>
	struct SeriesFromCfImpl
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedCf &cf, const ArgsTuple &argsTuple)
		{
			typedef typename Series::TermType::cf_type cf_series_type;
			cf_series_type cf_series;
			SeriesFromCfImpl<RequestedCf, typename cf_series_type::TermType::cf_type>::run(cf_series, cf,argsTuple);

			retval.insert(typename Series::TermType(cf_series, typename Series::TermType::key_type()), argsTuple);
		}
	};


	template <class Cf>
	struct SeriesFromCfImpl<Cf, Cf>
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Cf &cf, const ArgsTuple &argsTuple)
		{
			retval.insert(typename Series::TermType(cf, typename Series::TermType::key_type()), argsTuple);
		}
	};


	// This functor will disentangle and build the flattened terms iterating recursively through the echelon levels.
	template <int N>
	class SeriesFlattener
	{
		public:
			PIRANHA_STATIC_CHECK(N > 0,"");
			template <class CfSeries, class Term, class ArgsTuple>
			static void run(CfSeries &cf_series, Term &term, std::vector<Term> &out, const ArgsTuple &argsTuple)
			{
				// For each coefficient series (which is residing inside the insertion term), we create a copy of it,
				// then we insert one by one its original terms and, step by step, we go deeper into the recursion.
				PIRANHA_ASSERT(!cf_series.empty());
				const CfSeries tmp(cf_series);

				for (typename CfSeries::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
                {
					cf_series.clearTerms();
					cf_series.insert(*it, argsTuple);
					SeriesFlattener<N - 1>::run(cf_series.begin()->cf, term, out, argsTuple);
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

	template <class T, class Enable = void>
	struct BaseSeriesAddSelector
	{
		template <class Derived, class ArgsTuple>
		static Derived &run(Derived &series, const T &x, const ArgsTuple &argsTuple)
		{
			return series.template mergeWithNumber<true>(x, argsTuple);
		}
	};


	template <class T>
	struct BaseSeriesAddSelector<T, typename boost::enable_if<boost::is_base_of<BaseSeriesTag, T> >::type>
	{
		template <class Derived, class ArgsTuple>
		static Derived &run(Derived &series, const T &other, const ArgsTuple &argsTuple)
		{
			return series.template mergeTerms<true>(other, argsTuple);
		}
	};


	template <class T, class Enable = void>
	struct BaseSeriesSubtractSelector
	{
		template <class Derived, class ArgsTuple>
		static Derived &run(Derived &series, const T &x, const ArgsTuple &argsTuple)
		{
			return series.template mergeWithNumber<false>(x, argsTuple);
		}
	};


	template <class T>
	struct BaseSeriesSubtractSelector<T, typename boost::enable_if<boost::is_base_of<BaseSeriesTag, T> >::type>
	{
		template <class Derived, class ArgsTuple>
		static Derived &run(Derived &series, const T &other, const ArgsTuple &argsTuple)
		{
			return series.template mergeTerms<false>(other, argsTuple);
		}
	};


	template <class Derived, class T, class Enable = void>
	struct BaseSeriesMultiplySelector
	{
		template <class ArgsTuple>
		static Derived &run(Derived &series, const T &x, const ArgsTuple &argsTuple)
		{
			return series.multiplyCoefficientsBy(x, argsTuple);
		}
	};


	template <class Derived, class T>
	struct BaseSeriesMultiplySelector<Derived, T, typename boost::enable_if_c<boost::is_base_of<BaseSeriesTag, T>::value &&
		(boost::is_same<Derived, T>::value || boost::is_same<Derived, std::complex<T> >::value)>::type>
	{
		template <class ArgsTuple>
		static Derived &run(Derived &series, const T &other, const ArgsTuple &argsTuple)
		{
			series.multiply_by_series(other, argsTuple);
			return series;
		}
	};


	template <class T, class Enable = void>
	struct BaseSeriesEqualToSelector
	{
		template <class Derived>
		static bool run(const Derived &series, const T &x)
		{
			return series.genericNumericalComparison(x);
		}
	};


	template <class T>
	struct BaseSeriesEqualToSelector<T, typename boost::enable_if<boost::is_base_of<BaseSeriesTag, T> >::type>
	{
		template <class Derived>
		static bool run(const Derived &series, const T &other)
		{
			return series.genericSeriesComparison(other);
		}
	};
}

#endif
