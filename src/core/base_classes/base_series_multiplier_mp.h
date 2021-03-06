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

#ifndef PIRANHA_BASE_SERIES_MULTIPLIER_MP_H
#define PIRANHA_BASE_SERIES_MULTIPLIER_MP_H

#include <boost/tuple/tuple.hpp>

namespace piranha
{
	// Traverse the tuple of results of multiplication and insert each result into the multiplication result set.
	template <class ResultTuple>
	class InsertMultiplicationResult
	{
		public:
			template <class Series, class ArgsTuple>
			static void run(const ResultTuple &result, Series &series, const ArgsTuple &argsTuple)
			{
				series.insert(result.get_head(), argsTuple);
				InsertMultiplicationResult<typename ResultTuple::tail_type>::run(result.get_tail(), series, argsTuple);
			}
	};


	template <>
	class InsertMultiplicationResult<boost::tuples::null_type>
	{
		public:
			template <class Series, class ArgsTuple>
			static void run(const boost::tuples::null_type &, Series &, const ArgsTuple &) {}
	};
}

#endif
