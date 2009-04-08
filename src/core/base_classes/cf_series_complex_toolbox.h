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

#ifndef PIRANHA_CF_SERIES_COMPLEX_TOOLBOX_H
#define PIRANHA_CF_SERIES_COMPLEX_TOOLBOX_H

#include <complex>

#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	struct cf_series_complex {};

	template <>
	template <class RealDerived>
	class toolbox<cf_series_complex<RealDerived> >
	{
			typedef std::complex<RealDerived> Derived;
		public:
			template <class ArgsTuple>
			RealDerived real(const ArgsTuple &args_tuple) const {
				return derived_const_cast->base_real(args_tuple);
			}
			template <class ArgsTuple>
			RealDerived imag(const ArgsTuple &args_tuple) const {
				return derived_const_cast->base_imag(args_tuple);
			}
			template <class ArgsTuple>
			Derived &real(const RealDerived &r, const ArgsTuple &args_tuple) {
				return derived_cast->base_real(r,args_tuple);
			}
			template <class ArgsTuple>
			Derived &imag(const RealDerived &r, const ArgsTuple &args_tuple) {
				return derived_cast->base_imag(r,args_tuple);
			}
			template <class ArgsTuple>
			Derived &mult_by(const RealDerived &cs, const ArgsTuple &a) {
				return derived_cast->base_mult_by(cs,a);
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<double> &cx, const ArgsTuple &a) {
				return derived_cast->base_mult_by(cx,a);
			}
			template <class ArgsTuple>
			Derived &add(const RealDerived &cs, const ArgsTuple &a) {
				return derived_cast->base_add(cs,a);
			}
			template <class ArgsTuple>
			Derived &add(const std::complex<double> &cx, const ArgsTuple &a) {
				return derived_cast->base_add(cx,a);
			}
			template <class ArgsTuple>
			Derived &subtract(const RealDerived &cs, const ArgsTuple &a) {
				return derived_cast->base_subtract(cs,a);
			}
			template <class ArgsTuple>
			Derived &subtract(const std::complex<double> &cx, const ArgsTuple &a) {
				return derived_cast->base_subtract(cx,a);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<double> &cx, const ArgsTuple &a) {
				return derived_cast->base_divide_by(cx,a);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
