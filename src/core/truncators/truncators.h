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

#ifndef PIRANHA_TRUNCATORS_H
#define PIRANHA_TRUNCATORS_H

#include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/tuple/tuple.hpp>

#include "degree.h"
#include "norm.h"
#include "power_series.h"

namespace piranha
{
/// Namespace containing series truncators.
namespace truncators
{
	// Sequence of truncators type. NOTE: add here any new truncator.
	typedef boost::mpl::vector<norm,degree,power_series> truncator_types;

	// General functor for calling unset() on a truncator.
	template <class Truncator>
	struct truncator_unset {
		static void run()
		{
			Truncator::unset();
		}
	};

	// For power series truncators, there is not unset() function to be called.
	template <>
	struct truncator_unset<power_series> {
		static void run() {}
	};

	template <int N>
	struct unset_impl {
		static void run()
		{
			truncator_unset<typename boost::mpl::at<truncator_types,boost::mpl::int_<N> >::type>::run();
			unset_impl<N-1>::run();
		}
	};

	template <>
	struct unset_impl<-1> {
		static void run() {}
	};

	/// Unset all available truncators.
	inline void unset()
	{
		unset_impl<boost::mpl::size<truncator_types>::type::value - 1>::run();
	}
}
}

#endif
