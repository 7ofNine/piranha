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

#include <complex>
#include <iostream>
#include <string>

#include "../config.h"
#include "../math.h"
#include "../psym.h"
#include "../settings.h"
#include "../utils.h" // Lexical converter.
#include "numerical_container_complex_toolbox.h"
#include "series_builders.h"

// Convenience macros.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class T>
	class numerical_container_eval_type_determiner
	{
		public:
			typedef double type;
	};

	template <class T>
	class numerical_container_eval_type_determiner<std::complex<T> >
	{
		public:
			typedef std::complex<double> type;
	};

	/// Numerical container class.
	/**
	 * This class can be used as a base class for coefficients that consist of a
	 * numerical entity (double, GMP classes, etc.).
	 */
	template <class T, class Derived>
	class numerical_container
	{
			friend class numerical_container_complex_toolbox<Derived>;
		public:
			typedef typename numerical_container_eval_type_determiner<T>::type eval_type;
			typedef T numerical_type;
			template <class, class SubCachesCons, class>
			struct sub_cache_selector {
				typedef SubCachesCons type;
			};
			template <class, class SubCachesCons, class>
			struct ei_sub_cache_selector {
				typedef SubCachesCons type;
			};
			// Ctors.
			explicit numerical_container(): m_value(0) {}
			template <class ArgsTuple>
			explicit numerical_container(const std::string &s, const ArgsTuple &):
				m_value(utils::lexical_converter<T>(s)) {}
			template <class ArgsTuple>
			explicit numerical_container(const double &x, const ArgsTuple &): m_value(x) {}
			/// Ctor from psym.
			/**
			 * Sets m_value to one.
			 */
			template <class ArgsTuple>
			explicit numerical_container(const psym &, const int &, const ArgsTuple &): m_value(1) {}
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &) const {
				out_stream << m_value;
			}
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				print_plain(out_stream,args_tuple);
			}
			template <class ArgsTuple>
			void print_latex(std::ostream &out_stream, const ArgsTuple &) const {
// TODO: rework this.
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
			template <class TrimFlags>
			void trim_test(TrimFlags &) const {}
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &, const ArgsTuple &) const {
				return Derived(*derived_const_cast);
			}
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
			// If value is less than numerical zero in absolute value it is considered
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
			const T &eval(const double &, const ArgsTuple &) const {
				return m_value;
			}
			bool operator==(const Derived &other) const {
				return (m_value == other.m_value);
			}
			bool operator==(const double &x) const {
				return (m_value == x);
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
			Derived &mult_by(const double &x, const ArgsTuple &) {
				return mult_by_generic(x);
			}
			template <class ArgsTuple>
			Derived &mult_by(const Derived &x, const ArgsTuple &) {
				return mult_by_generic(x.value());
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
				return Derived(0, args_tuple);
			}
			template <class ArgsTuple>
			Derived inv(const ArgsTuple &args_tuple) const {
				return derived_const_cast->pow(-1,args_tuple);
			}
			template <class ArgsTuple>
			Derived root(const int &n, const ArgsTuple &args_tuple) const {
				if (n == 0) {
					throw division_by_zero();
				} else if (n == 1) {
					return Derived(*derived_const_cast);
				}
				return derived_const_cast->pow(1. / static_cast<double>(n), args_tuple);
			}
			template <class ArgsTuple>
			Derived besselJ(const int &, const ArgsTuple &) const {
				throw unsuitable("besselJ is not implemented for this coefficient type.");
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &, SubCaches &, const ArgsTuple &args_tuple) const {
				return numerical_cf_series_builder < boost::tuples::length<ArgsTuple>::value - 1 >::
					template run<RetSeries>(*derived_const_cast, args_tuple);
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &p, SubCaches &s, const ArgsTuple &a) const {
				return sub<RetSeries>(p,s,a);
			}
			/// Get value.
			const T &value() const {
				return m_value;
			}
			/// Set value.
			T &value() {
				return m_value;
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

#define NUMERICAL_CONTAINER_CTORS(derived) \
	explicit derived() {} \
	template <class ArgsTuple> \
	explicit derived(const std::string &s, const ArgsTuple &a): ancestor(s, a) {} \
	template <class ArgsTuple> \
	explicit derived(const double &val, const ArgsTuple &a): ancestor(val, a) {} \
	template <class ArgsTuple> \
	explicit derived(const piranha::psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}
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
