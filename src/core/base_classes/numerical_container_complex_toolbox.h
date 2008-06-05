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
#include "../integer_typedefs.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Toolbox for complex-specific methods of piranha::numerical_container.
	template <class realDerived>
	class numerical_container_complex_toolbox
	{
			typedef std::complex<realDerived> Derived;
			typedef realDerived value_type;
		public:
			// Ctors.
			numerical_container_complex_toolbox() {}
			explicit numerical_container_complex_toolbox(const std::complex<max_fast_int> &c) {
				derived_cast->m_value = c;
			}
			explicit numerical_container_complex_toolbox(const std::complex<double> &c) {
				derived_cast->m_value = c;
			}
			explicit numerical_container_complex_toolbox(const value_type &r) {
				derived_cast->m_value.real() = r.value();
			}
			explicit numerical_container_complex_toolbox(const value_type &r, const value_type &i) {
				derived_cast->m_value.real() = r.value();
				derived_cast->m_value.imag() = i.value();
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
			void real(const value_type &r, const ArgsTuple &) {
				derived_cast->m_value = r.value();
			}
			template <class ArgsTuple>
			void imag(const value_type &i, const ArgsTuple &) {
				derived_cast->m_value.imag() = i.value();
			}
			// Maths.
			template <class ArgsTuple>
			Derived &mult_by(const value_type &x, const ArgsTuple &) {
				return derived_cast->mult_by_generic(x.value());
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<max_fast_int> &c, const ArgsTuple &) {
				return derived_cast->mult_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<double> &c, const ArgsTuple &) {
				return derived_cast->mult_by_generic(c);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<max_fast_int> &c, const ArgsTuple &) {
				return complex_division_helper(c);
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<double> &c, const ArgsTuple &) {
				return complex_division_helper(c);
			}
			double display_data() const {
				return std::log10(1 + derived_const_cast->template norm<int>(0));
			}
		protected:
			template <class Number>
			Derived &complex_division_helper(const std::complex<Number> &x) {
				if (x.real() == 0 and x.imag() == 0) {
					throw division_by_zero();
				}
				// Do something only if we are not dividing by one.
				if (x.real() != 1 or x.imag() != 0) {
					Derived retval;
					const typename Derived::value_type::numerical_type abs2 = x.real() * x.real() + x.imag() * x.imag();
					retval.m_value.real() = derived_const_cast->m_value.real() * x.real() + derived_const_cast->m_value.imag() * x.imag();
					retval.m_value.imag() = derived_const_cast->m_value.imag() * x.real() - derived_const_cast->m_value.real() * x.imag();
					retval.m_value.real() /= abs2;
					retval.m_value.imag() /= abs2;
					derived_cast->m_value = retval.m_value;
				}
				return *derived_cast;
			}
	};

#define COMPLEX_NUMERICAL_CONTAINER_CTORS \
	template <class ArgsTuple> \
	explicit complex(const std::complex<piranha::max_fast_int> &c, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(c) {} \
	template <class ArgsTuple> \
	explicit complex(const std::complex<double> &c, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(c) {} \
	template <class ArgsTuple> \
	explicit complex(const value_type &r, const ArgsTuple &): complex_toolbox::numerical_container_complex_toolbox(r) {} \
	template <class ArgsTuple> \
	explicit complex(const value_type &r, const value_type &i, const ArgsTuple &): \
			complex_toolbox::numerical_container_complex_toolbox(r, i) {}
}

#undef derived_const_cast
#undef derived_cast

#endif
