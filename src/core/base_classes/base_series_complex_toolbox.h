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

#include "base_series.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class base_series_complex_toolbox
	{
			typedef std::complex<RealDerived> Derived;
		public:
			typedef RealDerived value_type;
			template <class ArgsTuple>
			RealDerived real(const ArgsTuple &args_tuple) const {
				return get_comp<0>(args_tuple);
			}
			template <class ArgsTuple>
			RealDerived imag(const ArgsTuple &args_tuple) const {
				return get_comp<1>(args_tuple);
			}
			template <class ArgsTuple>
			Derived &real(const RealDerived &r, const ArgsTuple &args_tuple) {
				subtract(real(args_tuple), args_tuple);
				add(r, args_tuple);
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &imag(const RealDerived &i, const ArgsTuple &args_tuple) {
				typedef typename RealDerived::const_iterator real_iterator;
				typedef typename Derived::term_type complex_term_type;
				complex_term_type tmp;
				// First let's remove the old imaginary part.
				RealDerived old_i(imag(args_tuple));
				const real_iterator old_i_it_f = old_i.end();
				for (real_iterator i_it = old_i.begin(); i_it != old_i_it_f; ++i_it) {
					tmp.m_key = i_it->m_key;
					tmp.m_cf.imag(i_it->m_cf, args_tuple);
					derived_cast->template insert<true, false>(tmp, args_tuple);
				}
				// Now add the new imaginary part.
				const real_iterator i_it_f = i.end();
				for (real_iterator i_it = i.begin(); i_it != i_it_f; ++i_it) {
					tmp.m_key = i_it->m_key;
					tmp.m_cf.imag(i_it->m_cf, args_tuple);
					derived_cast->insert(tmp, args_tuple);
				}
				return *derived_cast;
			}
			bool operator==(const std::complex<max_fast_int> &cn) const {
				return derived_const_cast->generic_numerical_comparison(cn);
			}
			bool operator==(const std::complex<double> &cx) const {
				return derived_const_cast->generic_numerical_comparison(cx);
			}
			bool operator!=(const std::complex<max_fast_int> &cn) const {
				return !(*this == cn);
			}
			bool operator!=(const std::complex<double> &cx) const {
				return !(*this == cx);
			}
			template <class ArgsTuple>
			Derived &add(const RealDerived &r, const ArgsTuple &args_tuple) {
				return derived_cast->template merge_terms<true>(r, args_tuple);
			}
			template <class ArgsTuple>
			Derived &subtract(const RealDerived &r, const ArgsTuple &args_tuple) {
				return derived_cast->template merge_terms<false>(r, args_tuple);
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<max_fast_int> &cn, const ArgsTuple &args_tuple) {
				return mult_by_complex(cn, args_tuple);
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<double> &cx, const ArgsTuple &args_tuple) {
				return mult_by_complex(cx, args_tuple);
			}
			template <class ArgsTuple>
			Derived &mult_by(const RealDerived &r, const ArgsTuple &args_tuple) {
				derived_cast->multiply_by_series(r, args_tuple);
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<max_fast_int> &cn, const ArgsTuple &args_tuple) {
				return divide_by_complex(cn, args_tuple);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<double> &cx, const ArgsTuple &args_tuple) {
				return divide_by_complex(cx, args_tuple);
			}
		protected:
			// Use N = 0 for real, N != 0 for imag.
			template <int N, class Real, class ArgsTuple>
			static Real get_cf_comp(const std::complex<Real> &c, const ArgsTuple &args_tuple) {
				if (N) {
					return c.imag(args_tuple);
				} else {
					return c.real(args_tuple);
				}
			}
			template <int N, class ArgsTuple>
			RealDerived get_comp(const ArgsTuple &args_tuple) const {
				typedef typename Derived::const_iterator complex_iterator;
				RealDerived retval;
				const complex_iterator c_it_f = derived_const_cast->end();
				for (complex_iterator c_it = derived_const_cast->begin(); c_it != c_it_f; ++c_it) {
					typename RealDerived::term_type tmp(get_cf_comp<N>(c_it->m_cf, args_tuple), c_it->m_key);
					retval.insert(tmp, args_tuple);
				}
				return retval;
			}
			template <class ArgsTuple>
			void construct_from_real(const RealDerived &r, const ArgsTuple &args_tuple) {
				// Make sure we are being called from an empty series.
				p_assert(derived_const_cast->empty());
				derived_cast->insert_range(r.begin(),r.end(),args_tuple);
			}
			template <class ArgsTuple>
			void construct_from_real_imag(const RealDerived &r, const RealDerived &i, const ArgsTuple &args_tuple) {
				typedef typename RealDerived::const_iterator real_iterator;
				typedef typename Derived::term_type complex_term_type;
				// Make sure we are being called from an empty series.
				p_assert(derived_const_cast->empty());
				// Let's build the real part first.
				construct_from_real(r, args_tuple);
				// Now let's proceed to the imaginary part.
				const real_iterator i_it_f = i.end();
				complex_term_type tmp;
				for (real_iterator i_it = i.begin(); i_it != i_it_f; ++i_it) {
					tmp.m_key = i_it->m_key;
					tmp.m_cf.imag(i_it->m_cf, args_tuple);
					derived_cast->insert(tmp, args_tuple);
				}
			}
		private:
			template <class Complex, class ArgsTuple>
			Derived &divide_by_complex(const Complex &c, const ArgsTuple &args_tuple) {
				if (c.real() == 0 && c.imag() == 0) {
					throw division_by_zero();
				} else if (c.real() != 1 || c.imag() != 0) {
					derived_cast->divide_coefficients_by(c, args_tuple);
				}
				return *derived_cast;
			}
			template <class Complex, class ArgsTuple>
			Derived &mult_by_complex(const Complex &c, const ArgsTuple &args_tuple) {
				if (c.real() == 0 && c.imag() == 0) {
					derived_cast->clear_terms();
				} else if (c.real() != 1 || c.imag() != 0) {
					derived_cast->multiply_coefficients_by(c, args_tuple);
				}
				return *derived_cast;
			}
	};

#define COMPLEX_E0_SERIES_TERM(term_name) term_name<std::complex<Cf>,Key,'|',Allocator>
#define COMPLEX_E0_SERIES(series_name) std::complex<E0_SERIES(series_name)>
#define COMPLEX_E0_SERIES_BASE_ANCESTOR(term_name,series_name) piranha::base_series<COMPLEX_E0_SERIES_TERM(term_name),'\n', \
	Allocator,COMPLEX_E0_SERIES(series_name) >

#define COMPLEX_E1_SERIES_TERM(term_name,cf_name) term_name<std::complex<cf_name>,Key1,'|',Allocator>
#define COMPLEX_E1_SERIES(series_name) std::complex<E1_SERIES(series_name)>
#define COMPLEX_E1_SERIES_BASE_ANCESTOR(term_name,cf_name,series_name) piranha::base_series<COMPLEX_E1_SERIES_TERM(term_name,cf_name), \
	'\n',Allocator,COMPLEX_E1_SERIES(series_name) >
}

#undef derived_const_cast
#undef derived_cast

#endif
