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

#ifndef PIRANHA_NAMED_SERIES_SPECIAL_FUNCTIONS_H
#define PIRANHA_NAMED_SERIES_SPECIAL_FUNCTIONS_H

#include "../integer_typedefs.h"
#include "base_series_special_functions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	struct named_series_special_functions: public base_series_special_functions<Derived> {
		Derived besselJ(const max_fast_int &order) const {
			Derived retval(derived_const_cast->b_besselJ(order, derived_const_cast->m_arguments));
			retval.m_arguments = derived_const_cast->m_arguments;
			return retval;
		}
		Derived dbesselJ(const max_fast_int &order) const {
			Derived retval(derived_const_cast->b_dbesselJ(order, derived_const_cast->m_arguments));
			retval.m_arguments = derived_const_cast->m_arguments;
			return retval;
		}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif