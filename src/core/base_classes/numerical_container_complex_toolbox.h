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

#ifndef PIRANHA_NUMERICAL_CONTAINER_COMPLEX_TOOLBOX_H
#define PIRANHA_NUMERICAL_CONTAINER_COMPLEX_TOOLBOX_H

#include <cmath>
#include <complex>

#include "../exceptions.h"
#include "../mp.h"
#include "../settings.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Toolbox for complex-specific methods of piranha::numerical_container.
	template <class RealDerived>
	class numerical_container_complex_toolbox
	{
			typedef std::complex<RealDerived> Derived;
			typedef RealDerived value_type;
		public:
			// NOTE: are these really good or it is better to employ the ctor macros below?
			numerical_container_complex_toolbox() {}
			explicit numerical_container_complex_toolbox(const std::complex<double> &c) {
				derived_cast->m_value = c;
			}
			explicit numerical_container_complex_toolbox(const std::complex<mp_rational> &c) {
				derived_cast->m_value = c;
			}
			explicit numerical_container_complex_toolbox(const std::complex<mp_integer> &c) {
				derived_cast->m_value = c;
			}
			explicit numerical_container_complex_toolbox(const value_type &r) {
				//derived_cast->m_value.real() = r.value();
				derived_cast->m_value = typename Derived::numerical_type(r.value());
			}
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &) const {
				out_stream << '(' << derived_const_cast->m_value.real();
				if (derived_const_cast->m_value.imag() >= 0) {
					out_stream << '+';
				}
				out_stream << derived_const_cast->m_value.imag() << "j)";
			}
			template <class ArgsTuple>
			void print_tex(std::ostream &out_stream, const ArgsTuple &) const {
				out_stream << "\\left(" << derived_const_cast->m_value.real();
				if (derived_const_cast->m_value.imag() >= 0) {
					out_stream << '+';
				}
				out_stream << derived_const_cast->m_value.imag() << "\\, \\imath\\right)";
			}
			// Getters and setters.
			template <class ArgsTuple>
			value_type real(const ArgsTuple &) const {
				value_type retval;
				retval.m_value = derived_const_cast->value().real();
				return retval;
			}
			template <class ArgsTuple>
			value_type imag(const ArgsTuple &) const {
				value_type retval;
				retval.m_value = derived_const_cast->value().imag();
				return retval;
			}
			template <class ArgsTuple>
			void set_real(const value_type &r, const ArgsTuple &) {
				derived_cast->m_value = r.value();
			}
			template <class ArgsTuple>
			void set_imag(const value_type &i, const ArgsTuple &) {
				// NOTE: this code works in gcc, but it is not standard.
				// derived_cast->m_value.imag() = i.value();
				derived_cast->m_value = typename Derived::numerical_type(derived_const_cast->m_value.real(),
					typename value_type::numerical_type(i.value()));
			}
			// Maths.
			template <class ArgsTuple>
			Derived &mult_by(const value_type &x, const ArgsTuple &) {
				return derived_cast->mult_by_generic(x.value());
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<double> &c, const ArgsTuple &) {
				return derived_cast->mult_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<mp_rational> &c, const ArgsTuple &) {
				return derived_cast->mult_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<mp_integer> &c, const ArgsTuple &) {
				return derived_cast->mult_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<double> &c, const ArgsTuple &) {
				return derived_cast->divide_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<mp_rational> &c, const ArgsTuple &) {
				return derived_cast->divide_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<mp_integer> &c, const ArgsTuple &) {
				return derived_cast->divide_by_generic(c);
			}
			// We must rewrite all comparisons, otherwise there will be conflicts with those defined in
			// numerical container.
			bool operator==(const Derived &x) const {
				return (derived_const_cast->m_value == x.m_value);
			}
			bool operator==(const RealDerived &x) const {
				return (derived_const_cast->m_value.real() == x.m_value && derived_const_cast->m_value.imag() == 0);
			}
			bool operator==(const double &x) const {
				return (derived_const_cast->m_value.real() == x && derived_const_cast->m_value.imag() == 0);
			}
			bool operator==(const std::complex<double> &cx) const {
				return (derived_const_cast->m_value == cx);
			}
			bool operator==(const mp_rational &q) const {
				return (derived_const_cast->m_value == q);
			}
			bool operator==(const std::complex<mp_rational> &cq) const {
				return (derived_const_cast->m_value == cq);
			}
			bool operator==(const mp_integer &z) const {
				return (derived_const_cast->m_value == z);
			}
			bool operator==(const std::complex<mp_integer> &cz) const {
				return (derived_const_cast->m_value == cz);
			}
	};

#define COMPLEX_NUMERICAL_CONTAINER_CTORS(...) \
	template <class ArgsTuple> \
	explicit complex(const std::complex<double> &c, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(c) {} \
	template <class ArgsTuple> \
	explicit complex(const std::complex<piranha::mp_rational> &qc, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(qc __VA_ARGS__) {} \
	template <class ArgsTuple> \
	explicit complex(const std::complex<piranha::mp_integer> &cz, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(cz __VA_ARGS__) {} \
	template <class ArgsTuple> \
	explicit complex(const value_type &r, const ArgsTuple &): complex_toolbox::numerical_container_complex_toolbox(r) {} \
	template <class ArgsTuple> \
	explicit complex(const value_type &r, const value_type &i, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(r, i) {} \
	using complex_toolbox::operator==; \
	using complex_toolbox::print_pretty; \
	using complex_toolbox::print_tex;
}

#undef derived_const_cast
#undef derived_cast

#endif
