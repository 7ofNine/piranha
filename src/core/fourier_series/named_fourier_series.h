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

#ifndef PIRANHA_NAMED_FOURIER_SERIES_H
#define PIRANHA_NAMED_FOURIER_SERIES_H

#include <string>
#include <utility>
#include <vector>

#include "../base_classes/toolbox.h"
#include "../exceptions.h"
#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	struct named_fourier_series {};

	/// Named Fourier series toolbox.
	template <class Derived>
	class toolbox<named_fourier_series<Derived> >
	{
		public:
			Derived integrate(const std::string &name) const
			{
				typedef typename ntuple<std::vector<std::pair<bool, size_t> >, 1>::type pos_tuple_type;
				const psym p(name);
				const pos_tuple_type pos_tuple = psyms2pos(vector_psym(1,p),derived_const_cast->m_arguments);
				Derived retval;
				if (pos_tuple.get_head()[0].first) {
					retval = derived_const_cast->base_integrate(pos_tuple,derived_const_cast->m_arguments);
					retval.m_arguments = derived_const_cast->m_arguments;
					retval.trim();
				} else {
					piranha_throw(value_error,"cannot integrate fourier series if integration variable is not present among the series' arguments");
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif