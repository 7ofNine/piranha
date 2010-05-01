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

#include <vector>

#include "../config.h"
#include "../exceptions.h"

namespace piranha
{
	// Implementation of factory methods from coefficients and keys.
	template <class RequestedKey, class SeriesKey>
	struct series_from_key_impl
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedKey &key, const ArgsTuple &args_tuple)
		{
			typedef typename Series::term_type::cf_type cf_series_type;
			cf_series_type cf_series;
			series_from_key_impl<RequestedKey,typename cf_series_type::term_type::key_type>::run(
				cf_series,key,args_tuple);
			retval.insert(typename Series::term_type(cf_series,typename Series::term_type::key_type()),args_tuple);
		}
	};

	template <class Key>
	struct series_from_key_impl<Key,Key>
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Key &key, const ArgsTuple &args_tuple)
		{
			retval.insert(typename Series::term_type(typename Series::term_type::cf_type(1,args_tuple),key),args_tuple);
		}
	};

	template <class RequestedCf, class SeriesCf>
	struct series_from_cf_impl
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedCf &cf, const ArgsTuple &args_tuple)
		{
			typedef typename Series::term_type::cf_type cf_series_type;
			cf_series_type cf_series;
			series_from_cf_impl<RequestedCf,typename cf_series_type::term_type::cf_type>::run(
				cf_series,cf,args_tuple);
			retval.insert(typename Series::term_type(cf_series,typename Series::term_type::key_type()),args_tuple);
		}
	};

	template <class Cf>
	struct series_from_cf_impl<Cf,Cf>
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Cf &cf, const ArgsTuple &args_tuple)
		{
			retval.insert(typename Series::term_type(cf,typename Series::term_type::key_type()),args_tuple);
		}
	};

	// This functor will disentangle and build the flattened terms iterating recursively through the echelon levels.
	template <int N>
	class series_flattener
	{
		public:
			p_static_check(N > 0,"");
			template <class CfSeries, class Term, class ArgsTuple>
			static void run(CfSeries &cf_series, Term &term, std::vector<Term> &out, const ArgsTuple &args_tuple)
			{
				// For each coefficient series (which is residing inside the insertion term), we create a copy of it,
				// then we insert one by one its original terms and, step by step, we go deeper into the recursion.
				piranha_assert(!cf_series.empty());
				const CfSeries tmp(cf_series);
				for (typename CfSeries::const_iterator it = tmp.begin(); it != tmp.end(); ++it) {
					cf_series.clear_terms();
					cf_series.insert(*it,args_tuple);
					series_flattener<N - 1>::run(cf_series.begin()->m_cf,term,out,args_tuple);
				}
			}
	};

	template <>
	class series_flattener<0>
	{
		public:
			template <class Cf, class Term, class ArgsTuple>
			static void run(Cf &, Term &term, std::vector<Term> &out, const ArgsTuple &)
			{
				out.push_back(term);
			}
	};
}

#endif
