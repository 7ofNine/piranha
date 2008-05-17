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

#ifndef PIRANHA_NUMERICAL_CONTAINER_H
#define PIRANHA_NUMERICAL_CONTAINER_H

#include <iostream>
#include <string>

#include "../arg_manager.h"
#include "../integer_typedefs.h"
#include "../psym.h"
#include "../utils.h" // Lexical converter.
#include "../type_traits.h"
#include "numerical_container_complex_toolbox.h"

// Convenience macros.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Numerical container class.
	/**
	 * This class can be used as a base class for coefficients that consist of a
	 * numerical entity (double, GMP classes, etc.).
	 */
	template <class T, class Derived>
	class numerical_container
	{
			// Alias for evaluation type.
			typedef typename eval_type<Derived>::type eval_type;
			friend class numerical_container_complex_toolbox<Derived>;
		public:
			// Start implementation of basic pseries coefficient interface.
			//------------
			// Ctors.
			explicit numerical_container(): m_value(0) {}
			template <class ArgsTuple>
			explicit numerical_container(const std::string &s, const ArgsTuple &): m_value(utils::lexical_converter<T>(s)) {}
			template <class ArgsTuple>
			explicit numerical_container(const max_fast_int &n, const ArgsTuple &): m_value(n) {}
			template <class ArgsTuple>
			explicit numerical_container(const double &x, const ArgsTuple &): m_value(x) {}
			/// Ctor from psym.
			/**
			 * Sets m_value to one.
			 */
			template <class ArgsTuple>
			explicit numerical_container(const psym_p &, const int &, const ArgsTuple &): m_value(1) {}
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &) const {
				stream_manager::setup_print(out_stream);
				out_stream << m_value;
			}
			template <class ArgsTuple>
			void print_latex(std::ostream &out_stream, const ArgsTuple &) const {
// TODO: rework this.
//         stream_manager::setup_print(out_stream);
//         out_stream << "$" << m_value << "$";
			}
			// Manipulation
			Derived &swap(Derived &dc) {
				std::swap(m_value, dc.m_value);
				return *derived_cast;
			}
			template <class ArgsTuple>
			void pad_right(const ArgsTuple &) {}
			template <class ArgsTuple, class Layout>
			void apply_layout(const ArgsTuple &, const Layout &) {}
			// Probing.
			size_t atoms() const {
				return 1;
			}
			template <class ArgsTuple>
			bool checkup(const ArgsTuple &) const {
				return true;
			}
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(m_value);
			}
			static const size_t max_size = 0;
			// If value is less than numericalzero  in absolute value it is considered
			// to be zero.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &a) const {
				return (static_cast<Derived const *>(this)->norm(a) < settings::numerical_zero());
			}
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &) const {
				return true;
			}
			template <class ArgsTuple>
			bool needs_padding(const ArgsTuple &) const {
				return false;
			}
			template <class ArgsTuple>
			const eval_type &eval(const double &, const ArgsTuple &) const {
				return m_value;
			}
			// Maths.
			template <class ArgsTuple>
			void invert_sign(const ArgsTuple &) {
				m_value *= -1;
			}
			template <class ArgsTuple>
			Derived &add(const Derived &val2, const ArgsTuple &) {
				return add_generic(val2.value());
			}
			template <class ArgsTuple>
			Derived &subtract(const Derived &val2, const ArgsTuple &) {
				return subtract_generic(val2.value());
			}
			template <class ArgsTuple>
			Derived &mult_by(const max_fast_int &n, const ArgsTuple &) {
				return mult_by_generic(n);
			}
			template <class ArgsTuple>
			Derived &mult_by(const double &x, const ArgsTuple &) {
				return mult_by_generic(x);
			}
			template <class ArgsTuple>
			Derived &mult_by(const Derived &x, const ArgsTuple &) {
				return mult_by_generic(x.value());
			}
			template <class ArgsTuple>
			Derived &divide_by(const max_fast_int &n, const ArgsTuple &) {
				return divide_by_generic(n);
			}
			template <class ArgsTuple>
			Derived &divide_by(const double &x, const ArgsTuple &) {
				return divide_by_generic(x);
			}
			// Multiply and add.
			template <class Derived2, class ArgsTuple>
			void addmul(const Derived &x1, const Derived2 &x2, const ArgsTuple &) {
				m_value += x1.m_value * x2.value();
			}
			template <class PosTuple, class ArgsTuple>
			Derived partial(const PosTuple &, const ArgsTuple &args_tuple) const {
				return Derived((max_fast_int)0, args_tuple);
			}
			/// Get value.
			const T &value() const {
				return m_value;
			}
			/// Set value.
			T &value() {
				return m_value;
			}
			max_fast_int degree() const {
				return 0;
			}
			max_fast_int min_degree() const {
				return 0;
			}
			template <class Vector>
			void upload_min_exponents(const Vector &) const {}
			template <class Vector>
			void test_min_exponents(const Vector &) const {}
			template <class Vector, class ArgsTuple>
			bool test_expo_limits(const Vector &, const ArgsTuple &) const {
				return true;
			}
		protected:
			template <class U>
			Derived &assign_self(const U &x) {
				m_value = x.value();
				return *derived_cast;
			}
			template <class U>
			Derived &add_generic(const U &x) {
				m_value += x;
				return *derived_cast;
			}
			template <class U>
			Derived &subtract_generic(const U &x) {
				m_value -= x;
				return *derived_cast;
			}
			template <class U>
			Derived &mult_by_generic(const U &x) {
				m_value *= x;
				return *derived_cast;
			}
			template <class U>
			Derived &divide_by_generic(const U &x) {
				m_value /= x;
				return *derived_cast;
			}
		protected:
			// Data member.
			T m_value;
	};
}

namespace std
{
	// Overloads for I/O operators.
	template <class T, class Derived>
	inline istream &operator>>(istream &is, piranha::numerical_container<T, Derived> &nc)
	{
		string tmp;
		getline(is, tmp);
		nc = piranha::numerical_container<T, Derived>(piranha::utils::lexical_converter<T>(tmp));
		return is;
	}

	template <class T, class Derived>
	inline ostream &operator<<(ostream &os, const piranha::numerical_container<T, Derived> &nc)
	{
		os << nc.value();
		return os;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
