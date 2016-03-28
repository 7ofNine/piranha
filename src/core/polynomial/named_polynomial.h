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

#include "../Psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Named polynomial toolbox.
	template <class Derived>
	class NamedPolynomial
	{
		public:

			Derived integrate(std::string const &name) const
			{
				typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, 1>::Type PositionTupleType;
				const Psym p(name);
                // get positon of symbol in the arguments tuple
				const PositionTupleType positionTuple = psyms2pos(VectorPsym(1, p), derived_const_cast->arguments());
				Derived retval;
				
				if (positionTuple.get_head()[0].first) 
				{
					retval = derived_const_cast->baseIntegrate(positionTuple, derived_const_cast->arguments());
					retval.setArguments(derived_const_cast->arguments());
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
