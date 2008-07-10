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

#ifndef PIRANHA_SERIES_BUILDERS_H
#define PIRANHA_SERIES_BUILDERS_H

#include <boost/tuple/tuple.hpp>

#include "../config.h"

namespace piranha {
	template <int ArgsIndex>
	class numerical_cf_series_builder
	{
			p_static_check(ArgsIndex > 0, "Invalid arguments tuple index in series builder from numerical coefficient.");
		public:
			template <class Series, class NumericalCf, class ArgsTuple>
			static Series run(const NumericalCf &ncf, const ArgsTuple &args_tuple) {
				typedef typename Series::term_type term_type;
				typedef typename term_type::cf_type cf_type;
				typedef typename term_type::key_type key_type;
				Series retval;
				retval.insert(
					term_type(numerical_cf_series_builder<ArgsIndex - 1>::template run<cf_type>(ncf,args_tuple),key_type()),
					args_tuple
				);
				return retval;
			}
	};

	template <>
	class numerical_cf_series_builder<0>
	{
		public:
			template <class Series, class NumericalCf, class ArgsTuple>
			static Series run(const NumericalCf &ncf, const ArgsTuple &args_tuple) {
				typedef typename Series::term_type term_type;
				typedef typename term_type::key_type key_type;
				Series retval;
				retval.insert(term_type(ncf,key_type()),args_tuple);
				return retval;
			}
	};

	template <int N, int Pos>
	class series_from_key_helper
	{
		public:
			template <class Series, class Key, class ArgsTuple>
			static Series run(const Key &key, const ArgsTuple &args_tuple) {
				typedef typename Series::term_type term_type;
				typedef typename term_type::cf_type cf_type;
				typedef typename term_type::key_type key_type;
				Series retval;
				retval.insert(
					term_type(series_from_key_helper<N-1,Pos>::template run<cf_type>(key,args_tuple),key_type()),
					args_tuple
				);
				return retval;
			}
	};

	template <int N>
	class series_from_key_helper<N,N>
	{
		public:
			template <class Series, class Key, class ArgsTuple>
			static Series run(const Key &key, const ArgsTuple &args_tuple) {
				typedef typename Series::term_type term_type;
				typedef typename term_type::cf_type cf_type;
				Series retval;
				retval.insert(
					term_type(cf_type(static_cast<max_fast_int>(1),args_tuple),key),
					args_tuple
				);
				return retval;
			}
	};

	class key_series_builder
	{
		public:
			template <class Series, class Key, class ArgsTuple>
			static Series run(const Key &key, const ArgsTuple &args_tuple) {
				p_static_check(boost::tuples::length<ArgsTuple>::value > 0, "Zero-size arguments tuple in builder from key.");
				p_static_check(boost::tuples::length<ArgsTuple>::value >= Key::position,
					"Key position in arguments tuple is greater than arguments' tuple size.");
				return series_from_key_helper<boost::tuples::length<ArgsTuple>::value-1,Key::position>::
					template run<Series>(key,args_tuple);
			}
	};
}

#endif
