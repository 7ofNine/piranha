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

#ifndef PIRANHA_NAMED_SERIES_MP_H
#define PIRANHA_NAMED_SERIES_MP_H

#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

#include "base_series_tag.h"

namespace piranha
{
	template <class T, class Enable = void>
	struct NamedSeriesAddSelector
	{
		template <class Derived>
		static Derived & run(Derived &series, const T &x)
		{
			return series.base_add(x, series.argumentsTuple);
		}
	};

	template <class T>
	struct NamedSeriesAddSelector<T, typename boost::enable_if<boost::is_base_of<BaseSeriesTag, T> >::type>
	{
		template <class Derived>
		static Derived & run(Derived & series, const T & other)
		{
			return series.template merge_with_series<true>(other);
		}
	};



	template <class T, class Enable = void>
	struct NamedSeriesSubtractSelector
	{
		template <class Derived>
		static Derived & run(Derived &series, const T &x)
		{
			return series.base_subtract(x, series.argumentsTuple);
		}
	};

	template <class T>
	struct NamedSeriesSubtractSelector<T,typename boost::enable_if<boost::is_base_of<BaseSeriesTag,T> >::type>
	{
		template <class Derived>
		static Derived &run(Derived &series, const T &other)
		{
			return series.template merge_with_series<false>(other);
		}
	};
	 

	// select which type of multiplicaiton to execute
	// multiply with a consant
	template <class T, class Enable = void>
	struct NamedSeriesMultiplySelector
	{
		template <class Derived>
		static Derived &run(Derived &series, const T &x)
		{
			return series.mult_number_helper(x);
		}
	};
	 
	//multiply with another series
	template <class T>
	struct NamedSeriesMultiplySelector<T, typename boost::enable_if<boost::is_base_of<BaseSeriesTag, T> >::type>
	{
		template <class Derived>
		static Derived & run(Derived & series, const T & other)
		{
			return series.mult_by_series(other);
		}
	};
	 
	template <class T, class Enable = void>
	struct NamedSeriesEqualitySelector
	{
		template <class Derived>
		static bool run(const Derived & series, const T & x)
		{
			return series.base_equal_to(x);
		}
	};

	template <class T>
	struct NamedSeriesEqualitySelector<T,typename boost::enable_if<boost::is_base_of<BaseSeriesTag, T> >::type>
	{
		template <class Derived>
		static bool run(const Derived &series, const T &other)
		{
			return series.series_comparison(other);
		}
	};
}

#endif
