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

#ifndef PIRANHA_NAMED_POWER_SERIES_H
#define PIRANHA_NAMED_POWER_SERIES_H

#include <vector>

#include "../ntuple.h"
#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	struct named_power_series {};

	/// Named power series toolbox.
	template <class Derived>
	class toolbox<named_power_series<Derived> >
	{
		public:
			int partial_degree(const std::vector<psym> &v_psym) const {
				typedef typename ntuple<std::vector<size_t>, Derived::n_arguments_sets>::type pos_tuple_type;
				typedef typename Derived::args_tuple_type args_tuple_type;
				pos_tuple_type pos_tuple;
				v_psym_2_pos_tuple<args_tuple_type>::run(derived_const_cast->m_arguments);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
