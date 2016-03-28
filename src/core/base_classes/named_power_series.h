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

#include <cstddef>
#include <string>
#include <vector>

#include "../Psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Named power series toolbox.
	template <class Degree, class Derived>
	class NamedPowerSeries
	{
		public:

			Degree partialDegree(std::vector<std::string> const &names) const
			{
				return derived_const_cast->base_partial_degree(psyms2pos(names2psyms(names), derived_const_cast->arguments()));
			}


			Degree partialOrder(const std::vector<std::string> &names) const
			{
                // that is what names2psyms does. 
				//VectorPsym symbols;
				//symbols.reserve(symbols.size());
				//for (std::size_t i = 0; i < vs.size(); ++i) 
                //{
				//	symbols.push_back(Psym(names[i]));
				//}

				return derived_const_cast->base_partial_order(psyms2pos(names2psyms(names), derived_const_cast->arguments()));
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
