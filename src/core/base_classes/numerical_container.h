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
#include "../mp.h"
#include "../psym.h"
#include "../settings.h"
#include "../utils.h" // Lexical converter.
#include "numerical_container_complex_toolbox.h"

// Convenience macros.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class T>
	struct numerical_container_eval_type_determiner
	{
		typedef double type;
	};

	template <class T>
	struct numerical_container_eval_type_determiner<std::complex<T> >
	{
		typedef std::complex<double> type;
	};

	template <class T, class U>
	struct in_place_transform {
		static const U &run(const U &u)
		{
			return u;
		}
	};

	template <class U>
	struct in_place_transform<double,U> {
		static double run(const U &u)
		{
			return u.to_double();
		}
	};

	template <>
	struct in_place_transform<double,double> {
		static const double &run(const double &u)
		{
			return u;
		}
	};

	template <class U>
	struct in_place_transform<std::complex<double>,U> {
		static std::complex<double> run(const U &u)
		{
			return u.to_double();
		}
	};

	template <>
	struct in_place_transform<std::complex<double>,std::complex<double> > {
		static const std::complex<double> &run(const std::complex<double> &u)
		{
			return u;
		}
	};

	template <>
	struct in_place_transform<std::complex<double>,double> {
		static const double &run(const double &u)
		{
			return u;
		}
	};

	/// Numerical container class.
	/**
	 * This class can be used as a base class for coefficients that consist of a
	 * numerical entity (double, MP classes, etc.).
	 */
	template <class T, class Derived>
	class numerical_container
	{
			friend class numerical_container_complex_toolbox<Derived>;
		public:
			/// Typedef for evaluation type.
			typedef typename numerical_container_eval_type_determiner<T>::type eval_type;
			/// Alias for internal type.
			typedef T numerical_type;
			template <class, class SubCachesCons, class>
			struct sub_cache_selector {
				typedef SubCachesCons type;
			};
			template <class, class SubCachesCons, class>
			struct ei_sub_cache_selector {
				typedef SubCachesCons type;
			};
			/// Default constructor, initialises internal value to 0.
			explicit numerical_container(): m_value(0) {}
			/// Constructor from string.
			/**
			 * Will call boost::lexical_converter internally.
			 */
			template <class ArgsTuple>
			explicit numerical_container(const std::string &s, const ArgsTuple &):
				m_value(utils::lexical_converter<T>(s)) {}
			/// Constructor from double.
			template <class ArgsTuple>
			explicit numerical_container(const double &x, const ArgsTuple &): m_value(x) {}
			/// Constructor from piranha::mp_rational.
			template <class ArgsTuple>
			explicit numerical_container(const mp_rational &q, const ArgsTuple &): m_value(q) {}
			/// Constructor from piranha::mp_integer.
			template <class ArgsTuple>
			explicit numerical_container(const mp_integer &z, const ArgsTuple &): m_value(z) {}
			/// Ctor from psym.
			/**
			 * Sets internal value to one.
			 */
			template <class ArgsTuple>
			explicit numerical_container(const psym &, const int &, const ArgsTuple &): m_value(1) {}
			/// Print in plain mode.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &) const {
				out_stream << m_value;
			}
			/// Print in pretty mode. Equivalent to print_plain.
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				print_plain(out_stream,args_tuple);
			}
			/// Swap content using std::swap.
			Derived &swap(Derived &dc) {
				std::swap(m_value, dc.m_value);
				return *derived_cast;
			}
			/// Pad right.
			template <class ArgsTuple>
			void pad_right(const ArgsTuple &) {}
			/// Apply layout.
			template <class Layout, class ArgsTuple>
			void apply_layout(const Layout &, const ArgsTuple &) {}
			/// Test if trimming is possible.
			template <class TrimFlags>
			void trim_test(TrimFlags &) const {}
			/// Trim.
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &, const ArgsTuple &) const {
				return Derived(*derived_const_cast);
			}
			/// Number of atoms. Returns 1.
			size_t atoms() const {
				return 1;
			}
			/// Norm. Returns std::abs of internal value.
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(m_value);
			}
			/// Test for ignorability.
			/**
			 * Returns true if norm() is less than settings::numerical_zero().
			 */
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &a) const {
				return (static_cast<Derived const *>(this)->norm(a) < settings::numerical_zero());
			}
			/// Insertability test. Returns true.
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &) const {
				return true;
			}
			/// Padding test. Returns false.
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
			bool operator==(const mp_rational &q) const {
				return (m_value == q);
			}
			bool operator==(const mp_integer &z) const {
				return (m_value == z);
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
			Derived &mult_by(const mp_rational &q, const ArgsTuple &) {
				return mult_by_generic(q);
			}
			template <class ArgsTuple>
			Derived &mult_by(const mp_integer &z, const ArgsTuple &) {
				return mult_by_generic(z);
			}
			template <class ArgsTuple>
			Derived &mult_by(const Derived &x, const ArgsTuple &) {
				return mult_by_generic(x.value());
			}
			template <class ArgsTuple>
			Derived &divide_by(const double &x, const ArgsTuple &) {
				return divide_by_generic(x);
			}
			template <class ArgsTuple>
			Derived &divide_by(const mp_rational &q, const ArgsTuple &) {
				return divide_by_generic(q);
			}
			template <class ArgsTuple>
			Derived &divide_by(const mp_integer &z, const ArgsTuple &) {
				return divide_by_generic(z);
			}
			// Multiply and add.
			template <class Derived2, class ArgsTuple>
			void addmul(const Derived &x1, const Derived2 &x2, const ArgsTuple &) {
				m_value += x1.m_value * x2.value();
			}
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &, const ArgsTuple &) const {
				return Series();
			}
			template <class ArgsTuple>
			Derived inv(const ArgsTuple &args_tuple) const {
				return derived_const_cast->pow(-1,args_tuple);
			}
			template <class ArgsTuple>
			Derived root(const int &n, const ArgsTuple &args_tuple) const {
				if (n == 0) {
					piranha_throw(zero_division_error,"cannot calculate zero-th root");
				} else if (n == 1) {
					return Derived(*derived_const_cast);
				}
				return derived_const_cast->pow(1. / static_cast<double>(n), args_tuple);
			}
			template <class ArgsTuple>
			Derived besselJ(const int &, const ArgsTuple &) const {
				piranha_throw(not_implemented_error,
					"besselJ is not implemented for this coefficient type");
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &, SubCaches &, const ArgsTuple &args_tuple) const {
				return RetSeries::base_series_from_cf(*derived_const_cast, args_tuple);
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
				m_value += in_place_transform<T,U>::run(x);
				return *derived_cast;
			}
			template <class U>
			Derived &subtract_generic(const U &x) {
				m_value -= in_place_transform<T,U>::run(x);
				return *derived_cast;
			}
			template <class U>
			Derived &mult_by_generic(const U &x) {
				m_value *= in_place_transform<T,U>::run(x);
				return *derived_cast;
			}
			template <class U>
			Derived &divide_by_generic(const U &x) {
				m_value /= in_place_transform<T,U>::run(x);
				return *derived_cast;
			}
		protected:
			T m_value;
	};

	#define NUMERICAL_CONTAINER_CTORS(derived,...) \
	explicit derived(): ancestor() {} \
	template <class ArgsTuple> \
	explicit derived(const std::string &s, const ArgsTuple &a): ancestor(s, a) {} \
	template <class ArgsTuple> \
	explicit derived(const double &val, const ArgsTuple &a): ancestor(val, a) {} \
	template <class ArgsTuple> \
	explicit derived(const piranha::mp_rational &val, const ArgsTuple &a): ancestor(val __VA_ARGS__ , a) {} \
	template <class ArgsTuple> \
	explicit derived(const piranha::mp_integer &val, const ArgsTuple &a): ancestor(val __VA_ARGS__ , a) {} \
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
