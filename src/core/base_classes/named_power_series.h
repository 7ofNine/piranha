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

#include "../psym.h"
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	struct named_power_series {};

	/// Named power series toolbox.
	template <class Derived>
	struct toolbox<named_power_series<Derived> >
	{
		int partial_degree(const vector_psym &v) const {
			return derived_const_cast->base_partial_degree(psyms2pos(v,derived_const_cast->m_arguments));
		}
		int partial_min_degree(const vector_psym &v) const {
			return derived_const_cast->base_partial_min_degree(psyms2pos(v,derived_const_cast->m_arguments));
		}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
