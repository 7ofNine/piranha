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

#include "../settings.h"
#include "named_series.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class named_series_complex
	{
			typedef std::complex<RealDerived> Derived;
		public:
			RealDerived real() const {
				RealDerived retval(derived_const_cast->base_real(derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}
			RealDerived imag() const {
				RealDerived retval(derived_const_cast->base_imag(derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}
			void set_real(const RealDerived &r) {
				derived_cast->merge_args(r);
				derived_cast->base_set_real(r, derived_const_cast->arguments());
				derived_cast->trim();
			}
			void set_imag(const RealDerived &i) {
				derived_cast->merge_args(i);
				derived_cast->base_set_imag(i, derived_const_cast->arguments());
				derived_cast->trim();
			}
			RealDerived abs() const {
				RealDerived retval = derived_const_cast->base_abs(derived_const_cast->arguments());
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}
			RealDerived abs2() const {
				RealDerived retval = derived_const_cast->base_abs2(derived_const_cast->arguments());
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}
			Derived conjugate() const {
				Derived retval = derived_const_cast->base_conjugate(derived_const_cast->arguments());
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


		protected:

			void construct_from_real(const RealDerived &r) {
				derived_cast->setArguments(r.arguments());
				derived_cast->base_construct_from_real(r, derived_cast->arguments());
				derived_cast->trim();
			}


			void construct_from_real_imag(const RealDerived &r, const RealDerived &i)
			{
				derived_cast->setArguments(r.arguments());
				derived_cast->merge_args(i);
				// NOTE: here we have to copy r and merge its arguments with i, otherwise
				// there will be arguments mismatches. All series, this, r and i must share the
				// same arguments layout.
				RealDerived r_copy(r);
				r_copy.merge_args(i);
				derived_cast->base_construct_from_real_imag(r_copy, i, derived_cast->arguments());
				derived_cast->trim();
			}
	};

#define COMPLEX_E0_SERIES_NAMED_ANCESTOR(args,term_name,series_name) piranha::NamedSeries<args,term_name,COMPLEX_E0_SERIES(series_name)>

#define COMPLEX_NAMED_SERIES_CTORS(real_series) \
	explicit complex(const complex<double> &cx) { \
		*this = base_series_from_number(cx,this->arguments()); \
		this->trim(); \
	} \
	explicit complex(const real_series &r) { \
		this->construct_from_real(r); \
	} \
	explicit complex(const real_series &r, const real_series &i) { \
		this->construct_from_real_imag(r, i); \
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
