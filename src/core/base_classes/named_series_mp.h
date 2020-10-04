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

#include "../type_traits.h"
#include "base_series_tag.h"


#include <type_traits>

namespace piranha
{
	template <class T>
	struct NamedSeriesAddSelector
	{

		template <class Derived>
		static Derived & run(Derived &series, T const &x)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.template mergeWithSeries<true>(x);
			}
			else
			{
				return series.baseAdd(x, series.argumentsTuple);
			}
		}
	};


	template <class T>
	struct NamedSeriesSubtractSelector
	{
		template <class Derived>
		static Derived& run(Derived& series, T const& x)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.template mergeWithSeries<false>(x);
			}
			else
			{
				return series.baseSubtract(x, series.argumentsTuple);
			}
		}
	};


	// select which type of multiplicaiton to execute
	// multiply with a consant
	template <typename T>
	struct NamedSeriesMultiplySelector
	{
		template <class Derived>
		static Derived& run(Derived& series, T const& x)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.multiplyBySeries(x);
			}
			else
			{
				return series.multiplyNumberHelper(x);
			}
		}
	};


	template <typename T>
	struct NamedSeriesEqualitySelector
	{
		template <class Derived>
		static bool run(const Derived& series, T const& x)
		{
			if constexpr (PiranhaSeries<T>)
			{
				return series.isEqualTo(x);
			}
			else
			{
				return series.baseEqualTo(x);
			}
		}
	};

}

#endif
