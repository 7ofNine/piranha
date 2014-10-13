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

#ifndef PIRANHA_NAMED_HARMONIC_SERIES_H
#define PIRANHA_NAMED_HARMONIC_SERIES_H

#include <string>
#include <vector>

#include "../Psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Named harmonic series toolbox.
	template <class HDegree, class Derived>
	class named_harmonic_series
	{
		public:

			HDegree partialHarmonicDegree(const std::vector<std::string> &vs) const
			{
				return derived_const_cast->basePartialHarmonicDegree(psyms2pos(names2psyms(vs), derived_const_cast->arguments()));
			}


			HDegree partialHarmonicOrder(const std::vector<std::string> &vs) const
			{
				VectorPsym v;
				v.reserve(vs.size());
				for (size_t i = 0; i < vs.size(); ++i) 
                {
					v.push_back(Psym(vs[i]));
				}

				return derived_const_cast->basePartialHarmonicOrder(psyms2pos(names2psyms(vs), derived_const_cast->arguments()));
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
