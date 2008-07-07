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

#ifndef PIRANHA_COMMON_FOURIER_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_FOURIER_SERIES_TOOLBOX_H

#include <complex>

#include "../poisson_series_common/jacobi_anger_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class common_fourier_series_toolbox: public jacobi_anger_toolbox<0, Derived>
	{
			typedef jacobi_anger_toolbox<0, Derived> jacang_ancestor;
		public:
			std::complex<Derived> complexp() const {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				typedef typename complex_term_type::key_type key_type;
				std::complex<Derived> retval;
				if (derived_const_cast->is_single_cf()) {
					retval.insert(complex_term_type(derived_const_cast->template nth_index<0>().begin()->
						m_cf.complexp(derived_const_cast->m_arguments),key_type()),
						derived_const_cast->m_arguments);
				} else {
					// Expand using Jacobi-Anger's identity.
					jacang_ancestor::jacobi_anger(derived_const_cast->template nth_index<0>().end(),
													retval, derived_const_cast->m_arguments);
					retval.m_arguments = derived_const_cast->m_arguments;
				}
				return retval;
			}
			Derived cos() const {
				return complexp().real();
			}
			Derived sin() const {
				return complexp().imag();
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
