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

#ifndef PIRANHA_NAMED_POLYNOMIAL_H
#define PIRANHA_NAMED_POLYNOMIAL_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Named polynomial toolbox.
	template <class Derived>
	class named_polynomial
	{
		public:
			Derived integrate(const std::string &name) const
			{
				typedef typename ntuple<std::vector<std::pair<bool, std::size_t> >, 1>::type pos_tuple_type;
				const psym p(name);
				const pos_tuple_type pos_tuple = psyms2pos(vector_psym(1, p), derived_const_cast->arguments());
				Derived retval;
				
				if (pos_tuple.get_head()[0].first) 
				{
					retval = derived_const_cast->base_integrate(pos_tuple,derived_const_cast->arguments());
					retval.set_arguments(derived_const_cast->arguments());
					retval.trim();
				} else 
				{
					retval = *derived_const_cast;
					retval *= Derived(p);
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
