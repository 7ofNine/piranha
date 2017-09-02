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

#ifndef PIRANHA_NTUPLE_H
#define PIRANHA_NTUPLE_H

#include <boost/tuple/tuple.hpp>

#include "config.h"

namespace piranha
{
	// Wrapper for tuple of homogeneous types.
	//used for argument types with T= std::vector<Psym>
	template <class T, int N>
	struct NTuple
	{
        static_assert(N > 0, "");

		typedef boost::tuples::cons < T, typename NTuple < T, N - 1 >::Type > Type;
	};

	template <class T>
	struct NTuple<T, 1> 
	{
		typedef boost::tuples::cons<T, boost::tuples::null_type> Type;
	};

	template <class T>
	struct NTuple<T, 0> 
	{
		typedef boost::tuples::null_type Type;
	};
}

#endif
