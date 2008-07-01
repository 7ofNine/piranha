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

#ifndef PIRANHA_SHARED_ARGS_MP_H
#define PIRANHA_SHARED_ARGS_MP_H

#include <boost/tuple/tuple.hpp>

namespace piranha
{
	template <class ArgsTuple, class StaticArgsTuple>
	struct shared_args_assign_tuple {
		static void run(const ArgsTuple &args_tuple, StaticArgsTuple &t) {
			t.get_head() = args_tuple.get_head();
			shared_args_assign_tuple<typename ArgsTuple::tail_type, typename StaticArgsTuple::tail_type>::
			run(args_tuple.get_tail(), t.get_tail());
		}
	};

	template <class StaticArgsTuple>
	struct shared_args_assign_tuple<boost::tuples::null_type, StaticArgsTuple> {
		static void run(const boost::tuples::null_type &, StaticArgsTuple &) {}
	};
}

#endif
