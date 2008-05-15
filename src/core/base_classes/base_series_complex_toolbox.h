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

#ifndef PIRANHA_BASE_SERIES_COMPLEX_TOOLBOX_H
#define PIRANHA_BASE_SERIES_COMPLEX_TOOLBOX_H

#include <complex>

#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class base_series_complex_toolbox
	{
			typedef std::complex<RealDerived> Derived;
		public:
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<max_fast_int> &cn, const ArgsTuple &args_tuple) {
				if (cn.real() == 0 and cn.imag() == 0) {
					Derived retval;
					derived_cast->swap_terms(retval);
				} else if (cn.real() != 1 or cn.imag() != 0) {
					Derived retval(derived_cast->multiply_coefficients_by(cn, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<double> &cx, const ArgsTuple &args_tuple) {
				if (cx.real() == 0 and cx.imag() == 0) {
					Derived retval;
					derived_cast->swap_terms(retval);
				} else if (cx.real() = ! 1 or cx.imag() != 0) {
					Derived retval(derived_cast->multiply_coefficients_by(cx, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &mult_by(const RealDerived &r, const ArgsTuple &args_tuple)
			{
				Derived retval(derived_cast->multiply_by_series(r, args_tuple));
				// Grab the terms accumulated into return value.
				derived_cast->swap_terms(retval);
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<max_fast_int> &cn, const ArgsTuple &args_tuple) {
				if (cn.real() == 0 and cn.imag() == 0) {
					throw division_by_zero();
				} else if (cn.real() != 1 or cn.imag() != 0) {
					Derived retval(derived_cast->divide_coefficients_by(cn, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<double> &cx, const ArgsTuple &args_tuple) {
				if (cx.real() == 0 and cx.imag() == 0) {
					throw division_by_zero();
				} else if (cx.real() != 1 or cx.imag() != 0) {
					Derived retval(derived_cast->divide_coefficients_by(cx, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
