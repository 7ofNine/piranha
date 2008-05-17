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

#include "base_series.h"
#include "base_series_complex_toolbox.h"
#include "named_series.h"
#include "../integer_typedefs.h"
#include "../settings.h"

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
			RealDerived real() const {
				RealDerived retval(ancestor::real(derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				return retval;
			}
			RealDerived imag() const {
				RealDerived retval(ancestor::imag(derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				return retval;
			}
			Derived &real(const RealDerived &r) {
				derived_cast->merge_args(r);
				return ancestor::real(r, derived_const_cast->m_arguments);
			}
			Derived &imag(const RealDerived &i) {
				derived_cast->merge_args(i);
				return ancestor::imag(i, derived_const_cast->m_arguments);
			}
			Derived &operator+=(const std::complex<max_fast_int> &cn) {
				return derived_cast->template merge_with_number<true>(cn, derived_cast->m_arguments);
			}
			Derived &operator+=(const std::complex<double> &cx) {
				return derived_cast->template merge_with_number<true>(cx, derived_cast->m_arguments);
			}
			Derived &operator+=(const RealDerived &r) {
				return derived_cast->template merge_with_series<true>(r);
			}
			Derived &operator-=(const std::complex<max_fast_int> &cn) {
				return derived_cast->template merge_with_number<false>(cn, derived_cast->m_arguments);
			}
			Derived &operator-=(const std::complex<double> &cx) {
				return derived_cast->template merge_with_number<false>(cx, derived_cast->m_arguments);
			}
			Derived &operator-=(const RealDerived &r) {
				return derived_cast->template merge_with_series<false>(r);
			}
			Derived &operator*=(const std::complex<max_fast_int> &cn) {
				return ancestor::mult_by(cn, derived_const_cast->m_arguments);
			}
			Derived &operator*=(const std::complex<double> &cx) {
				return ancestor::mult_by(cx, derived_const_cast->m_arguments);
			}
			Derived &operator*=(const RealDerived &r) {
				return derived_cast->mult_by_series(r);
			}
			Derived &operator/=(const std::complex<max_fast_int> &cn) {
				return ancestor::divide_by(cn, derived_const_cast->m_arguments);
			}
			Derived &operator/=(const std::complex<double> &cx) {
				return ancestor::divide_by(cx, derived_const_cast->m_arguments);
			}
		protected:
			void construct_from_real(const RealDerived &r) {
				derived_cast->m_arguments = r.m_arguments;
				ancestor::construct_from_real(r, derived_cast->m_arguments);
			}
			void construct_from_real_imag(const RealDerived &r, const RealDerived &i) {
				derived_cast->m_arguments = r.m_arguments;
				derived_cast->merge_args(i);
				ancestor::construct_from_real_imag(r, i, derived_cast->m_arguments);
			}
	};

#define COMPLEX_NAMED_SERIES_TERM(term_name) term_name<std::complex<Cf>,Key,'|',Allocator>
#define COMPLEX_NAMED_SERIES(series_name) std::complex<REAL_NAMED_SERIES(series_name)>
#define COMPLEX_NAMED_SERIES_BASE_ANCESTOR(term_name,series_name) piranha::base_series<COMPLEX_NAMED_SERIES_TERM(term_name),'\n', \
	Allocator,COMPLEX_NAMED_SERIES(series_name) >
#define COMPLEX_NAMED_SERIES_NAMED_ANCESTOR(args,series_name) piranha::named_series<args,COMPLEX_NAMED_SERIES(series_name) >

#define COMPLEX_NAMED_SERIES_CTORS(complex_toolbox) \
	explicit complex(const complex<piranha::max_fast_int> &cn) { \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		base_ancestor::construct_from_number(cn,named_ancestor::m_arguments); \
	} \
	explicit complex(const complex<double> &cx) { \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		base_ancestor::construct_from_number(cx,named_ancestor::m_arguments); \
	} \
	explicit complex(const value_type &r) { \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		complex_toolbox::construct_from_real(r); \
	} \
	explicit complex(const value_type &r, const value_type &i) { \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		complex_toolbox::construct_from_real_imag(r, i); \
	}

}

#undef derived_const_cast
#undef derived_cast

#endif
