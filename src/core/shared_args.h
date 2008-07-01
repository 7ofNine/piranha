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

#ifndef PIRANHA_SHARED_ARGS_H
#define PIRANHA_SHARED_ARGS_H

#include "shared_args_mp.h"
#include "config.h"
#include "ntuple.h"
#include "psym.h"

namespace piranha
{
	// TODO: update doc.
	/// Manager for arguments.
	/**
	 * This class is used to manage information about arguments in those context where such information
	 * is not available. For instance when managing the elements of the multiindex container in a series
	 * class the functors used to order the terms do not know anything about arguments, since arguments do not
	 * appear inside each term. In those cases, before manipulating a multiindex container, an
	 * arg_manager::arg_assigner instance should be created, so that proper arguments are made available
	 * through the get() method.
	 */
	class __PIRANHA_VISIBLE shared_args
	{
			typedef ntuple<vector_psym_p, 10>::type temp_tuple_type;
		public:
			template <class ArgsTuple>
			static void set(const ArgsTuple &args_tuple) {
				shared_args_assign_tuple<ArgsTuple, temp_tuple_type>::run(args_tuple, m_args_tuple);
			}
			/// Retrieve const reference to arguments tuple.
			static const temp_tuple_type &get() {
				return m_args_tuple;
			}
		private:
			shared_args() {}
			~shared_args() {}
		private:
			static temp_tuple_type      m_args_tuple;
	};
}

#endif
