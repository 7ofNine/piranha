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

#ifndef PIRANHA_NUMERICAL_CONTAINER_COMPLEX_TOOLBOX_H
#define PIRANHA_NUMERICAL_CONTAINER_COMPLEX_TOOLBOX_H

#include <cmath>
#include <complex>

#include "../exceptions.h"
#include "../mp.h"
#include "../settings.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Toolbox for complex-specific methods of piranha::numerical_container.
	template <class RealDerived>
	class numerical_container_complex_toolbox
	{
		public:
			typedef std::complex<RealDerived> Derived;
			typedef RealDerived value_type;
			// Getters and setters.
			template <class ArgsTuple>
			value_type real(const ArgsTuple &) const {
				value_type retval;
				retval.set_value(derived_const_cast->get_value().real());
				return retval;
			}
			template <class ArgsTuple>
			value_type imag(const ArgsTuple &) const {
				value_type retval;
				retval.set_value(derived_const_cast->get_value().imag());
				return retval;
			}
			template <class ArgsTuple>
			void set_real(const value_type &r, const ArgsTuple &) {
				derived_cast->set_value(r.get_value());
			}
			template <class ArgsTuple>
			void set_imag(const value_type &i, const ArgsTuple &) {
				// NOTE: this code works in gcc, but it is not standard.
				// derived_cast->m_value.imag() = i.value();
				derived_cast->set_value(typename Derived::numerical_type(derived_const_cast->get_value().real(),
					typename value_type::numerical_type(i.get_value())));
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
