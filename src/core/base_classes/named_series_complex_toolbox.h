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
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

// TODO: here we are still missing all interoperability with complex mps.

namespace piranha
{
	template <class RealDerived>
	struct named_series_complex {};

	template <class RealDerived>
	class toolbox<named_series_complex<RealDerived> >
	{
			typedef std::complex<RealDerived> Derived;
		public:
			RealDerived real() const {
				RealDerived retval(derived_const_cast->base_real(derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			RealDerived imag() const {
				RealDerived retval(derived_const_cast->base_imag(derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			void set_real(const RealDerived &r) {
				derived_cast->merge_args(r);
				derived_cast->base_set_real(r, derived_const_cast->m_arguments);
				derived_cast->trim();
			}
			void set_imag(const RealDerived &i) {
				derived_cast->merge_args(i);
				derived_cast->base_set_imag(i, derived_const_cast->m_arguments);
				derived_cast->trim();
			}
			RealDerived abs() const {
				RealDerived retval = derived_const_cast->base_abs(derived_const_cast->m_arguments);
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			RealDerived abs2() const {
				RealDerived retval = derived_const_cast->base_abs2(derived_const_cast->m_arguments);
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			Derived conjugate() const {
				Derived retval = derived_const_cast->base_conjugate(derived_const_cast->m_arguments);
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			// NOTE: this is probably implementable in base fashion, by doing all the arguments merging dance
			// and then comparing. But since we do only one comparison (real() == x) in which all the dancing is already
			// done, we would not gain anything.
			bool operator==(const RealDerived &x) const
			{
				return (derived_const_cast->base_imag(derived_const_cast->m_arguments).empty() && real() == x);
			}
			// TODO: remove != operators when using boost::operators.
			bool operator!=(const RealDerived &x) const
			{
				return !(operator==(x));
			}
			bool operator==(const std::complex<double> &cx) const {
				return derived_const_cast->base_equal_to_complex_number(cx);
			}
			bool operator!=(const std::complex<double> &cx) const {
				return !(*this == cx);
			}
			Derived &operator+=(const std::complex<double> &cx) {
				return derived_cast->template merge_number_helper<true>(cx);
			}
			Derived &operator+=(const RealDerived &r) {
				return derived_cast->template merge_with_series<true>(r);
			}
			Derived &operator-=(const std::complex<double> &cx) {
				return derived_cast->template merge_number_helper<false>(cx);
			}
			Derived &operator-=(const RealDerived &r) {
				return derived_cast->template merge_with_series<false>(r);
			}
			Derived &operator*=(const std::complex<double> &cx) {
				return derived_cast->mult_number_helper(cx);
			}
			Derived &operator*=(const RealDerived &r) {
				return derived_cast->mult_by_series(r);
			}
			Derived &operator/=(const std::complex<double> &cx) {
				return derived_cast->divide_number_helper(cx);
			}
		protected:
			void construct_from_real(const RealDerived &r) {
				derived_cast->m_arguments = r.m_arguments;
				derived_cast->base_construct_from_real(r, derived_cast->m_arguments);
				derived_cast->trim();
			}
			void construct_from_real_imag(const RealDerived &r, const RealDerived &i) {
				derived_cast->m_arguments = r.m_arguments;
				derived_cast->merge_args(i);
				derived_cast->base_construct_from_real_imag(r, i, derived_cast->m_arguments);
				derived_cast->trim();
			}
	};

#define COMPLEX_E0_SERIES_NAMED_ANCESTOR(args,term_name,series_name) piranha::toolbox<piranha::named_series<args,term_name,COMPLEX_E0_SERIES(series_name) > >

#define COMPLEX_NAMED_SERIES_CTORS(real_series) \
	explicit complex(const complex<double> &cx) { \
		*this = base_series_from_number(cx,this->m_arguments); \
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
