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

#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../poisson_series_common/jacobi_anger_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class ArgsTuple>
	class term_norm_comparison
	{
		public:
			term_norm_comparison(const ArgsTuple &args_tuple):m_args_tuple(args_tuple) {}
			template <class Term>
			bool operator()(const Term &t1, const Term &t2) const {
				return t1.m_cf.norm(m_args_tuple) > t2.m_cf.norm(m_args_tuple);
			}
		private:
			const ArgsTuple &m_args_tuple;
	};

	template <class Derived>
	class common_fourier_series_toolbox:
		public jacobi_anger_toolbox<0, Derived>,
		public binomial_exponentiation_toolbox<Derived,term_norm_comparison>
	{
			typedef jacobi_anger_toolbox<0, Derived> jacang_ancestor;
		public:
			std::complex<Derived> complexp() const {
				std::complex<Derived> retval(complexp(derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			template <class ArgsTuple>
			std::complex<Derived> complexp(const ArgsTuple &args_tuple) const {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				typedef typename complex_term_type::key_type key_type;
				typedef typename Derived::const_iterator::type const_iterator;
				std::complex<Derived> retval;
				if (derived_const_cast->is_single_cf()) {
					retval.insert(complex_term_type(derived_const_cast->begin()->
													m_cf.complexp(args_tuple), key_type()),
								  args_tuple);
				} else {
					// Let's find out if there is a constant term.
					const_iterator it = derived_const_cast->begin();
					const const_iterator it_f = derived_const_cast->end();
					for (; it != it_f; ++it) {
						if (it->m_key.is_unity()) {
							break;
						}
					}
					// Expand using Jacobi-Anger's identity.
					jacang_ancestor::jacobi_anger(it, retval, args_tuple);
					if (it != it_f) {
						std::complex<Derived> tmp;
						tmp.insert(complex_term_type(it->m_cf.complexp(args_tuple),key_type()),args_tuple);
						retval.mult_by(tmp,args_tuple);
					}
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
