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

#include "../base_classes/toolbox.h"
#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class HDegree, class Derived>
	struct named_harmonic_series {};

	/// Named harmonic series toolbox.
	template <class HDegree, class Derived>
	struct toolbox<named_harmonic_series<HDegree,Derived> >
	{
		HDegree partial_h_degree(const std::vector<std::string> &vs) const
		{
			return derived_const_cast->base_partial_h_degree(psyms2pos(names2psyms(vs),derived_const_cast->m_arguments));
		}
		HDegree partial_h_order(const std::vector<std::string> &vs) const
		{
			vector_psym v;
			v.reserve(vs.size());
			for (size_t i = 0; i < vs.size(); ++i) {
				v.push_back(psym(vs[i]));
			}
			return derived_const_cast->base_partial_h_order(psyms2pos(names2psyms(vs),derived_const_cast->m_arguments));
		}
		/// Flip flavour.
		Derived flip_flavour() const
		{
			Derived retval;
			derived_const_cast->base_flip_flavour(retval,derived_const_cast->m_arguments);
			retval.m_arguments = derived_const_cast->m_arguments;
			retval.trim();
			return retval;
		}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif