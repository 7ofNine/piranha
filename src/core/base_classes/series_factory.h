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

#ifndef PIRANHA_SERIES_FACTORY_H
#define PIRANHA_SERIES_FACTORY_H

#include "toolbox.h"

namespace piranha
{
	template <class RequestedKey, class SeriesKey>
	struct series_from_key_impl_tag {};

	template <class RequestedKey, class SeriesKey>
	struct toolbox<series_from_key_impl_tag<RequestedKey,SeriesKey> >
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedKey &key, const ArgsTuple &args_tuple)
		{
			typedef typename Series::term_type::cf_type cf_series_type;
			cf_series_type cf_series;
			toolbox<series_from_key_impl_tag<RequestedKey,typename cf_series_type::term_type::key_type> >::run(
				cf_series,key,args_tuple);
			retval.insert(typename Series::term_type(cf_series,typename Series::term_type::key_type()),args_tuple);
		}
	};

	template <class Key>
	struct toolbox<series_from_key_impl_tag<Key,Key> >
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Key &key, const ArgsTuple &args_tuple)
		{
			retval.insert(typename Series::term_type(typename Series::term_type::cf_type(1,args_tuple),key),args_tuple);
		}
	};

	template <class RequestedKey, class SeriesKey>
	struct series_from_cf_impl_tag {};

	template <class RequestedCf, class SeriesCf>
	struct toolbox<series_from_cf_impl_tag<RequestedCf,SeriesCf> >
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const RequestedCf &cf, const ArgsTuple &args_tuple)
		{
			typedef typename Series::term_type::cf_type cf_series_type;
			cf_series_type cf_series;
			toolbox<series_from_cf_impl_tag<RequestedCf,typename cf_series_type::term_type::cf_type> >::run(
				cf_series,cf,args_tuple);
			retval.insert(typename Series::term_type(cf_series,typename Series::term_type::key_type()),args_tuple);
		}
	};

	template <class Cf>
	struct toolbox<series_from_cf_impl_tag<Cf,Cf> >
	{
		template <class Series, class ArgsTuple>
		static void run(Series &retval, const Cf &cf, const ArgsTuple &args_tuple)
		{
			retval.insert(typename Series::term_type(cf,typename Series::term_type::key_type()),args_tuple);
		}
	};
}

#endif
