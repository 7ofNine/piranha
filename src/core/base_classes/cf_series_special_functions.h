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

#ifndef PIRANHA_CF_SERIES_SPECIAL_FUNCTIONS_H
#define PIRANHA_CF_SERIES_SPECIAL_FUNCTIONS_H

#include <string>

#include "../exceptions.h"
#include "../math.h"
#include "../p_assert.h"
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	struct cf_series_special_functions {};

	template <>
	template <class Derived>
	class toolbox<cf_series_special_functions<Derived> >
	{
		public:
			/// Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived besselJ(const int &order, const ArgsTuple &args_tuple) const {
				return derived_const_cast->base_besselJ(order,args_tuple);
			}
			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived dbesselJ(const int &order, const ArgsTuple &args_tuple) const {
				return derived_const_cast->base_dbesselJ(order,args_tuple);
			}
			/// Bessel function of the first kind of integer order divided by its argument**m.
			template <class ArgsTuple>
			Derived besselJ_div_m(const int &order, const int &m, const ArgsTuple &args_tuple) const {
				return derived_const_cast->base_besselJ_div_m(order,m,args_tuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
