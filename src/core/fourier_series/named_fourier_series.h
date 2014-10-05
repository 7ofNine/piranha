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

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "../exceptions.h"
#include "../Psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Named Fourier series toolbox.
	template <class Derived>
	class named_fourier_series
	{
		public:
			Derived integrate(const std::string &name) const
			{
				typedef typename Ntuple<std::vector<std::pair<bool, std::size_t> >, 1>::type pos_tuple_type;
				const Psym p(name);
				const pos_tuple_type pos_tuple = psyms2pos(VectorPsym(1,p),derived_const_cast->arguments());
				Derived retval;
				if (pos_tuple.get_head()[0].first) {
					retval = derived_const_cast->base_integrate(pos_tuple,derived_const_cast->arguments());
					retval.set_arguments(derived_const_cast->arguments());
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
