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

#ifndef PIRANHA_NAMED_SERIES_COMPLEX_TOOLBOX_H
#define PIRANHA_NAMED_SERIES_COMPLEX_TOOLBOX_H

#include <complex>

#include "base_series_complex_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class named_series_complex_toolbox : public base_series_complex_toolbox<RealDerived>
	{
			typedef base_series_complex_toolbox<RealDerived> ancestor;
			typedef std::complex<RealDerived> Derived;
		public:
			Derived &operator*=(const std::complex<max_fast_int> &cn) {
				return ancestor::mult_by(cn, derived_const_cast->m_arguments);
			}
			Derived &operator*=(const std::complex<double> &cx) {
				return ancestor::mult_by(cx, derived_const_cast->m_arguments);
			}
			Derived &operator/=(const std::complex<max_fast_int> &cn) {
				return ancestor::divide_by(cn, derived_const_cast->m_arguments);
			}
			Derived &operator/=(const std::complex<double> &cx) {
				return ancestor::divide_by(cx, derived_const_cast->m_arguments);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
