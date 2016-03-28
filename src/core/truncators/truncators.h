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
	typedef boost::mpl::vector<Norm, Degree, PowerSeries> TruncatorTypes;


	// General functor for calling unset() on a truncator.
	template <class Truncator>
	class TruncatorUnset
    {
        public: 
		
        static void run()
		{
			Truncator::unset();
		}
	};


	// For power series truncators, there is not unset() function to be called.
	template <>
	struct TruncatorUnset<PowerSeries>
    {
		static void run() {}
	};


    // implementaton for unsetting all types of truncators.
    // i.e. truncation is set to inactive or level 0. This is truncator specific.
	template <int N>
	class UnsetImpl
     {
        public:

		static void run()
		{
			TruncatorUnset<typename boost::mpl::at<TruncatorTypes, boost::mpl::int_<N> >::type>::run();
			UnsetImpl<N-1>::run();
		}
	};

    // terminate recursion through truncator types
	template <>
	class UnsetImpl<-1>
    {
        public:

		static void run() {}
	};


	/// Unset all available truncators.
	inline void unset()
	{
		UnsetImpl<boost::mpl::size<TruncatorTypes>::type::value - 1>::run();
	}
}
}

#endif
