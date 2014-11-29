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

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class cf_series_complex
	{
			typedef std::complex<RealDerived> Derived;
		public:
			template <class ArgsTuple>
			RealDerived real(const ArgsTuple &argsTuple) const {
				return derived_const_cast->baseReal(argsTuple);
			}
			template <class ArgsTuple>
			RealDerived imag(const ArgsTuple &argsTuple) const {
				return derived_const_cast->baseImag(argsTuple);
			}
			template <class ArgsTuple>
			void set_real(const RealDerived &r, const ArgsTuple &argsTuple) {
				derived_cast->baseSetReal(r,argsTuple);
			}
			template <class ArgsTuple>
			void set_imag(const RealDerived &r, const ArgsTuple &argsTuple) {
				derived_cast->baseSetImag(r,argsTuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
