/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redis\bute it and/or modify  *
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

#ifndef PIRANHA_PIRANHA_GMP_H
#define PIRANHA_PIRANHA_GMP_H

#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/operators.hpp>
#include <cmath>
#include <complex>
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../utils.h"
#include "complex_generic_mp_container.h"

namespace piranha
{
	// Forward declaration of classes.
	class mp_rational;
	class mp_integer;

	#define FORWARDING_CTOR_DECL(class_type,arg_type) \
	/** \brief Constructor from arg_type. */ \
	explicit class_type(const arg_type &);
	#define FORWARDING_CTOR(class_type,arg_type) \
	inline class_type::class_type(const arg_type &x):m_value(x.get_internal()) {}
	#define FORWARDING_MATH_OPERATOR_DECL(class_type,arg_type,op) \
	/** \brief Operator op against arg_type. */ \
	class_type &operator op(const arg_type &);
	#define FORWARDING_MATH_OPERATOR(class_type,arg_type,op,...) \
	inline class_type & class_type::operator op(const arg_type &other) \
	{ \
		m_value op other __VA_ARGS__; \
		return *this; \
	}
	#define FORWARDING_DIVISION_OPERATOR_DECL(class_type,arg_type) \
	/** \brief Operator /= against arg_type. */ \
	class_type & operator /=(const arg_type &);
	#define FORWARDING_DIVISION_OPERATOR(class_type,arg_type,...) \
	inline class_type & class_type::operator /=(const arg_type &other) \
	{ \
		if (other == 0) { \
			piranha_throw(zero_division_error,"cannot divide by zero"); \
		} \
		m_value /= other __VA_ARGS__; \
		return *this; \
	}
	#define FORWARDING_COMPARISON_OPERATOR_DECL(class_type,arg_type,op) \
	/** \brief Operator op against arg_type. */ \
	bool operator op(const arg_type &) const;
	#define FORWARDING_COMPARISON_OPERATOR(class_type,arg_type,op,...) \
	inline bool class_type::operator op(const arg_type &other) const \
	{ \
		return (m_value op other __VA_ARGS__ ); \
	}

	/// Multiprecision rational class.
	/**
	 * Wraps a GMP mpq_class.
	 */
	class mp_rational:
		boost::ordered_field_operators<mp_rational,
		boost::ordered_field_operators<mp_rational, int,
		boost::ordered_field_operators<mp_rational, double,
		boost::ordered_field_operators<mp_rational, mp_integer
		> > > >
	{
		public:
			/// Default constructor.
			/**
			 * Initialises value to zero.
			 */
			explicit mp_rational():m_value(0) {}
			/// Constructor from std::string.
			/**
			 * Will raise a value_error exception if string is not valid. Valid strings include
			 * integer numbers and rationals, even with negative signs (e.g., "3", "-10", "3/4", "5/-6", etc.).
			 * @throws value_error if string is invalid.
			 * @throws zero_division_error if string is valid but denominator is zero.
			 */
			explicit mp_rational(const std::string &str):m_value(0)
			{
				construct_from_string(str.c_str());
			}
			/// Constructor from C string.
			/**
			 * @see mp_rational(const std::string &).
			 */
			explicit mp_rational(const char *str):m_value(0)
			{
				construct_from_string(str);
			}
			/// Constructor from integer.
			explicit mp_rational(const int &n):m_value(n) {}
			/// Constructor from integer numerator and denominator.
			/**
			 * @throws zero_division_error if denominator is zero.
			 */
			explicit mp_rational(const int &n, const int &d):m_value(0)
			{
				// Guard against division by zero.
				if (d == 0) {
					piranha_throw(zero_division_error,"cannot create rational with zero as denominator");
				}
				mpq_class tmp(n,d);
				mpq_canonicalize(tmp.get_mpq_t());
				// Swap content (more efficient than copying).
				mpz_swap(mpq_numref(m_value.get_mpq_t()),mpq_numref(tmp.get_mpq_t()));
				mpz_swap(mpq_denref(m_value.get_mpq_t()),mpq_denref(tmp.get_mpq_t()));
			}
			/// Constructor from double.
			explicit mp_rational(const double &x):m_value(x) {}
			/// Const reference to internal type.
			const mpq_class &get_internal() const
			{
				return m_value;
			}
			// Interoperability with self and plain old types.
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_rational,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_rational,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_rational,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_rational,mp_rational)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_rational,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_rational,<)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,int,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,int,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,int,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,int,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_rational,int)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,int,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,int,<)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,int,>)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,double,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,double,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,double,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,double,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_rational,double)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,double,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,double,<)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,double,>)
			// Interoperability with mp_integer.
			FORWARDING_CTOR_DECL(mp_rational,mp_integer)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_rational,mp_integer)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_integer,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_integer,<)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_integer,>)
			/// Cast to double.
			/**
			 * Uses the get_d() GMP routine internally.
			 */
			double to_double() const
			{
				// NOTE: the natural way to do it would seem to be this, but reading the GMP docs
				// it seems like this could fail in horrible ways in case of overflows. Need to check.
				return m_value.get_d();
			}
			/// Convert to integer.
			/**
			 * @throws value_error if the denominator is not unitary.
			 * @throws std::overflow_error if the numerator overflows int type.
			 */
			int to_int() const
			{
				if (m_value.get_den() != 1) {
					piranha_throw(value_error,"cannot convert rational to integer if denominator is non-unitary");
				}
				if (!m_value.get_num().fits_slong_p()) {
					piranha_throw(std::overflow_error,"numerator is too large while converting rational to long integer");
				}
				return m_value.get_num().get_si();
			}
			/// Swap content.
			/**
			 * Internally uses the mpz_swap GMP function on numerator and denominator.
			 */
			void swap(mp_rational &other)
			{
				mpz_swap(mpq_numref(m_value.get_mpq_t()),mpq_numref(other.m_value.get_mpq_t()));
				mpz_swap(mpq_denref(m_value.get_mpq_t()),mpq_denref(other.m_value.get_mpq_t()));
			}
			/// Hash value.
			/**
			 * Internally uses boost::hash_combine on the GMP limbs of numerator and denominator.
			 */
			size_t hash() const
			{
				size_t retval = 0;
				const __mpz_struct *num = mpq_numref(m_value.get_mpq_t()), *den = mpq_denref(m_value.get_mpq_t());
				const size_t num_limb_size = std::abs(num->_mp_size), den_limb_size = std::abs(den->_mp_size);
				for (size_t i = 0; i < num_limb_size; ++i) {
					boost::hash_combine(retval, num->_mp_d[i]);
				}
				for (size_t i = 0; i < den_limb_size; ++i) {
					boost::hash_combine(retval, den->_mp_d[i]);
				}
				return retval;
			}
			/// Negate in-place.
			/**
			 * Implemented through the mpq_neg GMP function.
			 */
			void negate()
			{
				mpq_neg(m_value.get_mpq_t(),m_value.get_mpq_t());
			}
			/// Negated copy.
			/**
			 * Call negate() on a copy of this and return it.
			 * @see negate().
			 */
			mp_rational operator-() const
			{
				mp_rational retval(*this);
				retval.negate();
				return retval;
			}
			/// Exponentiation.
			/**
			 * @throws zero_division_error if power is negative and value is zero.
			 * @throws value_error if power cannot be calculated exactly.
			 */
			mp_rational pow(const double &y) const
			{
				if (utils::is_integer(y)) {
					return pow_int((int)y);
				} else {
				}	return pow_double(y);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n_ is zero.
			 * @throws value_error if root cannot be calculated exactly.
			 */
			mp_rational root(const int &n_) const
			{
				mp_rational retval;
				if (n_ == 0) {
					piranha_throw(zero_division_error,"cannot calculate zero-th root of rational number");
				} else if (n_ == 1) {
					retval = *this;
					return retval;
				}
				const size_t n = (n_ > 0) ? n_ : -n_;
				if (!mpz_root(mpq_numref(retval.m_value.get_mpq_t()),mpq_numref(m_value.get_mpq_t()),n) ||
					!mpz_root(mpq_denref(retval.m_value.get_mpq_t()),mpq_denref(m_value.get_mpq_t()),n)) {
					piranha_throw(value_error,"rational number is not an exact nth root");
				}
				// Better to canonicalise, for peace of mind.
				retval.m_value.canonicalize();
				if (n_ < 0) {
					// Let's guard against division by zero below.
					if (retval == 0) {
						piranha_throw(zero_division_error,"cannot calculate negative root of zero");
					}
					mpq_inv(retval.m_value.get_mpq_t(), mpq_class(retval.m_value).get_mpq_t());
				}
				return retval;
			}
		private:
			// Will throw value_error if string is invalid, zero_division_error if string is valid
			// but contains zero as denominator.
			void construct_from_string(const char *str)
			{
				std::string tmp(str);
				// First let's trim the input string.
				boost::trim(tmp);
				// Let's try to separate numerator and denominator.
				std::vector<std::string> split_v;
				boost::split(split_v,tmp,boost::is_any_of("/"));
				switch (split_v.size()) {
					case 1:
						// If there is no "/" character, we assume that string represent an integer number,
						// hence we try to build mpq_class from it.
						try {
							m_value = mpq_class(split_v[0]);
						} catch (const std::invalid_argument &) {
							piranha_throw(value_error,"invalid string input");
						}
						break;
					case 2:
						{
						mpz_class num(0), den(0);
						try {
							num = mpz_class(split_v[0]);
							den = mpz_class(split_v[1]);
						} catch (const std::invalid_argument &) {
							piranha_throw(value_error,"invalid string input");
						}
						if (den == 0) {
							piranha_throw(zero_division_error,"cannot create rational with zero as denominator");
						}
						mpq_class tmp_q(num,den);
						mpq_canonicalize(tmp_q.get_mpq_t());
						m_value = tmp_q;
						}
						break;
					default:
						piranha_throw(value_error,"invalid string input");
				}
			}
			mp_rational pow_int(const int &n) const {
				mp_rational retval;
				if (m_value == 0 && n < 0) {
					piranha_throw(zero_division_error,"cannot raise zero to negative power");
				}
				if (n < 0) {
					mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (unsigned long int)(-n));
					mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (unsigned long int)(-n));
					// We need to canonicalize, since negative numbers may have gone to the denominator.
					retval.m_value.canonicalize();
				} else {
					mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (unsigned long int)n);
					mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (unsigned long int)n);
				}
				return retval;
			}
			mp_rational pow_double(const double &y) const {
				mp_rational retval;
				// If negative, only 1^-something is reasonable.
				if (y < 0) {
					if (m_value == 0) {
						piranha_throw(zero_division_error,"cannot raise zero to negative power");
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						piranha_throw(value_error,"cannot raise rational number different from unity to "
							"negative real power");
					}
				} else if (y == 0) {
					// If y == 0, then x**0 == 1 for every x.
					retval.m_value = 1;
				} else {
					// If y > 0, we can accept only 0^y and 1^y.
					if (m_value == 0) {
						retval.m_value = 0;
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						piranha_throw(value_error,"cannot raise rational number different from unity to "
							"positive real power");
					}
				}
				return retval;
			}
		private:
			mpq_class m_value;
	};

	FORWARDING_MATH_OPERATOR(mp_rational,mp_rational,+=,.m_value)
	FORWARDING_MATH_OPERATOR(mp_rational,mp_rational,-=,.m_value)
	FORWARDING_MATH_OPERATOR(mp_rational,mp_rational,*=,.m_value)
	FORWARDING_DIVISION_OPERATOR(mp_rational,mp_rational,.m_value)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,mp_rational,==,.m_value)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,mp_rational,<,.m_value)
	FORWARDING_MATH_OPERATOR(mp_rational,int,=)
	FORWARDING_MATH_OPERATOR(mp_rational,int,+=)
	FORWARDING_MATH_OPERATOR(mp_rational,int,-=)
	FORWARDING_MATH_OPERATOR(mp_rational,int,*=)
	FORWARDING_DIVISION_OPERATOR(mp_rational,int)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,int,==)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,int,<)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,int,>)
	FORWARDING_MATH_OPERATOR(mp_rational,double,=)
	FORWARDING_MATH_OPERATOR(mp_rational,double,+=)
	FORWARDING_MATH_OPERATOR(mp_rational,double,-=)
	FORWARDING_MATH_OPERATOR(mp_rational,double,*=)
	FORWARDING_DIVISION_OPERATOR(mp_rational,double)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,double,==)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,double,<)
	FORWARDING_COMPARISON_OPERATOR(mp_rational,double,>)

	/// Overload of out stream operator<< for piranha::mp_rational.
	inline std::ostream &operator<<(std::ostream &o, const mp_rational &q)
	{
		o << q.get_internal();
		return o;
	}

	/// Overload in stream operator>> for piranha::mp_rational.
	inline std::istream &operator>>(std::istream &i, mp_rational &q)
	{
		std::string tmp_str;
		std::getline(i,tmp_str);
		mp_rational tmp(tmp_str);
		q.swap(tmp);
		return i;
	}

	/// Overload hash_value function for piranha::mp_rational.
	inline size_t hash_value(const mp_rational &q)
	{
                return q.hash();
	}

	/// Multiprecision integer class.
	/**
	 * Wraps a GMP mpz_class.
	 */
	class mp_integer:
		boost::ordered_field_operators<mp_integer,
		boost::ordered_field_operators<mp_integer, int,
		boost::ordered_field_operators<mp_integer, double
		> > >
	{
		public:
			/// Default constructor.
			/**
			 * Initialises value to zero.
			 */
			explicit mp_integer():m_value(0) {}
			/// Constructor from std::string.
			explicit mp_integer(const std::string &str):m_value(str) {}
			/// Constructor from C string.
			explicit mp_integer(const char *str):m_value(str) {}
			/// Constructor from integer.
			explicit mp_integer(const int &n):m_value(n) {}
			/// Constructor from double.
			explicit mp_integer(const double &x):m_value(x) {}
			/// Const reference to internal type.
			const mpz_class &get_internal() const
			{
				return m_value;
			}
			// Interoperability with plain numerical types.
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_integer,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_integer,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_integer,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_integer,mp_integer)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,mp_integer,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,mp_integer,<)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,int,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,int,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,int,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,int,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_integer,int)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,int,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,int,<)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,int,>)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,double,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,double,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,double,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,double,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_integer,double)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,double,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,double,<)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_integer,double,>)
			// Interoperability with mp_rational.
			FORWARDING_CTOR_DECL(mp_integer,mp_rational)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_integer,mp_rational)
			/// Cast to double.
			/**
			 * Uses the get_d() GMP routine internally.
			 */
			double to_double() const
			{
				// NOTE: the natural way to do it would seem to be this, but reading the GMP docs
				// it seems like this could fail in horrible ways in case of overflows. Need to check.
				return m_value.get_d();
			}
			/// Convert to integer.
			/**
			 * @throws std::overflow_error if the numerator overflows int type.
			 */
			int to_int() const
			{
				if (!m_value.fits_sint_p()) {
					piranha_throw(std::overflow_error,"multiprecision integer too big to be converted to integer");
				}
				return m_value.get_si();
			}
			/// Swap content.
			/**
			 * Internally uses the mpz_swap GMP function.
			 */
			void swap(mp_integer &other)
			{
				mpz_swap(m_value.get_mpz_t(),other.m_value.get_mpz_t());
			}
			/// Hash value.
			/**
			 * Internally uses boost::hash_combine on the GMP limbs of the number.
			 */
			size_t hash() const
			{
				size_t retval = 0;
				const __mpz_struct *ptr = m_value.get_mpz_t();
				const size_t limb_size = std::abs(ptr->_mp_size);
				for (size_t i = 0; i < limb_size; ++i) {
					boost::hash_combine(retval, ptr->_mp_d[i]);
				}
				return retval;
			}
			/// Negate in-place.
			/**
			 * Implemented through the mpz_neg GMP function.
			 */
			void negate()
			{
				mpz_neg(m_value.get_mpz_t(),m_value.get_mpz_t());
			}
			/// Negated copy.
			/**
			 * Call negate() on a copy of this and return it.
			 * @see negate().
			 */
			mp_integer operator-() const
			{
				mp_integer retval(*this);
				retval.negate();
				return retval;
			}
			/// Multiply and add.
			/**
			 * *this += z1 * z2.
			 */
			void addmul(const mp_integer &z1, const mp_integer &z2)
			{
				mpz_addmul(m_value.get_mpz_t(), z1.m_value.get_mpz_t(), z2.m_value.get_mpz_t());
			}
			/// Exponentiation.
			/**
			 * @throws zero_division_error if power is negative and value is zero.
			 * @throws value_error if power cannot be calculated exactly.
			 */
			mp_integer pow(const double &y) const
			{
				if (utils::is_integer(y)) {
					return pow_int((int)y);
				} else {
				}	return pow_double(y);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n_ is zero.
			 * @throws value_error if root cannot be calculated exactly.
			 */
			mp_integer root(const int &n_) const
			{
				mp_integer retval;
				if (n_ == 0) {
					piranha_throw(zero_division_error,"cannot calculate zero-th root");
				} else if (n_ == 1) {
					retval = *this;
					return retval;
				} else if (n_ < 0) {
					piranha_throw(value_error,"integer numbers different from unity cannot be arguments of negative root");
				}
				const size_t n = static_cast<size_t>(n_);
				if (!mpz_root(retval.m_value.get_mpz_t(),m_value.get_mpz_t(),n)) {
					piranha_throw(value_error,"integer coefficient is not an exact nth root");
				}
				return retval;
			}
		private:
			mp_integer pow_int(const int &n) const
			{
				mp_integer retval;
				// If negative, only 1^-something is reasonable.
				if (n < 0) {
					if (m_value == 0) {
						piranha_throw(zero_division_error,"cannot divide by zero");
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						piranha_throw(value_error,"cannot raise integer number different from unity to "
							"negative integer power");
					}
				} else {
					mpz_pow_ui(retval.m_value.get_mpz_t(), m_value.get_mpz_t(), (size_t)n);
				}
				return retval;
			}
			mp_integer pow_double(const double &y) const
			{
				mp_integer retval;
				// If negative, only 1^-something is reasonable.
				if (y < 0) {
					if (m_value == 0) {
						piranha_throw(zero_division_error,"cannot divide by zero");
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						piranha_throw(value_error,"cannot raise integer number different from unity to negative real power");
					}
					// If y == 0, then x**0 == 1 for every x.
				} else if (y == 0) {
					retval.m_value = 1;
					// If y > 0, we can accept only 0^y and 1^y.
				} else {
					if (m_value == 0) {
						retval.m_value = 0;
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						piranha_throw(value_error,"cannot raise integer number different from unity to positive real power");
					}
				}
				return retval;
			}
		private:
			mpz_class m_value;
	};

	FORWARDING_MATH_OPERATOR(mp_integer,mp_integer,+=,.m_value)
	FORWARDING_MATH_OPERATOR(mp_integer,mp_integer,-=,.m_value)
	FORWARDING_MATH_OPERATOR(mp_integer,mp_integer,*=,.m_value)
	FORWARDING_DIVISION_OPERATOR(mp_integer,mp_integer,.m_value)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,mp_integer,==,.m_value)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,mp_integer,<,.m_value)
	FORWARDING_MATH_OPERATOR(mp_integer,int,=)
	FORWARDING_MATH_OPERATOR(mp_integer,int,+=)
	FORWARDING_MATH_OPERATOR(mp_integer,int,-=)
	FORWARDING_MATH_OPERATOR(mp_integer,int,*=)
	FORWARDING_DIVISION_OPERATOR(mp_integer,int)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,int,==)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,int,<)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,int,>)
	FORWARDING_MATH_OPERATOR(mp_integer,double,=)
	FORWARDING_MATH_OPERATOR(mp_integer,double,+=)
	FORWARDING_MATH_OPERATOR(mp_integer,double,-=)
	FORWARDING_MATH_OPERATOR(mp_integer,double,*=)
	FORWARDING_DIVISION_OPERATOR(mp_integer,double)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,double,==)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,double,<)
	FORWARDING_COMPARISON_OPERATOR(mp_integer,double,>)

	/// Overload of out stream operator<< for piranha::mp_integer.
	inline std::ostream &operator<<(std::ostream &o, const mp_integer &q)
	{
		o << q.get_internal();
		return o;
	}

	/// Overload in stream operator>> for piranha::mp_integer.
	inline std::istream &operator>>(std::istream &i, mp_integer &z)
	{
		std::string tmp_str;
		std::getline(i,tmp_str);
		mp_integer tmp(tmp_str);
		z.swap(tmp);
		return i;
	}

	/// Overload hash_value function for piranha::mp_integer.
	inline size_t hash_value(const mp_integer &z)
	{
                return z.hash();
	}

	// Mixed operations between mp types.
	FORWARDING_CTOR(mp_rational,mp_integer)
	FORWARDING_MATH_OPERATOR(mp_rational,mp_integer,=, .get_internal())
	FORWARDING_MATH_OPERATOR(mp_rational,mp_integer,+=, .get_internal())
	FORWARDING_MATH_OPERATOR(mp_rational,mp_integer,-=, .get_internal())
	FORWARDING_MATH_OPERATOR(mp_rational,mp_integer,*=, .get_internal())
	FORWARDING_DIVISION_OPERATOR(mp_rational,mp_integer, .get_internal())
	FORWARDING_COMPARISON_OPERATOR(mp_rational,mp_integer,==, .get_internal())
	FORWARDING_COMPARISON_OPERATOR(mp_rational,mp_integer,<, .get_internal())
	FORWARDING_COMPARISON_OPERATOR(mp_rational,mp_integer,>, .get_internal())
	FORWARDING_CTOR(mp_integer,mp_rational)
	FORWARDING_MATH_OPERATOR(mp_integer,mp_rational,=, .get_internal())
	FORWARDING_MATH_OPERATOR(mp_integer,mp_rational,+=, .get_internal())
	FORWARDING_MATH_OPERATOR(mp_integer,mp_rational,-=, .get_internal())
	FORWARDING_MATH_OPERATOR(mp_integer,mp_rational,*=, .get_internal())
	FORWARDING_DIVISION_OPERATOR(mp_integer,mp_rational, .get_internal())

	#undef FORWARDING_CTOR_DECL
	#undef FORWARDING_CTOR
	#undef FORWARDING_MATH_OPERATOR_DECL
	#undef FORWARDING_MATH_OPERATOR
	#undef FORWARDING_DIVISION_OPERATOR_DECL
	#undef FORWARDING_DIVISION_OPERATOR
	#undef FORWARDING_COMPARISON_OPERATOR_DECL
	#undef FORWARDING_COMPARISON_OPERATOR
}

namespace std
{
	/// Overload standard swap function for piranha::mp_rational.
	/**
	 * Will use piranha::mp_rational::swap() internally.
	 * @see piranha::mp_rational::swap().
	 */
	inline void swap(piranha::mp_rational &q1, piranha::mp_rational &q2)
	{
		q1.swap(q2);
	}

	/// Overload standard power function for piranha::mp_rational and double argument.
	/**
	 * @see piranha::mp_rational::pow.
	 */
	inline piranha::mp_rational pow(const piranha::mp_rational &q, const double &y)
	{
		return q.pow(y);
	}

	/// Overload standard power function for piranha::mp_rational and int argument.
	/**
	 * @see piranha::mp_rational::pow.
	 */
	inline piranha::mp_rational pow(const piranha::mp_rational &q, const int &y)
	{
		return q.pow(y);
	}

	/// Overload standard swap function for piranha::mp_integer.
	/**
	 * Will use piranha::mp_integer::swap() internally.
	 * @see piranha::mp_integer::swap().
	 */
	inline void swap(piranha::mp_integer &z1, piranha::mp_integer &z2)
	{
		z1.swap(z2);
	}

	/// Overload standard power function for piranha::mp_integer and double argument.
	/**
	 * @see piranha::mp_integer::pow.
	 */
	inline piranha::mp_integer pow(const piranha::mp_integer &z, const double &y)
	{
		return z.pow(y);
	}

	/// Overload standard power function for piranha::mp_integer and int argument.
	/**
	 * @see piranha::mp_integer::pow.
	 */
	inline piranha::mp_integer pow(const piranha::mp_integer &z, const int &y)
	{
		return z.pow(y);
	}

	#define CTOR_DECL(class_type,arg_type) \
	/** Constructor from arg_type. */ \
	explicit class_type(const arg_type &);
	#define CTOR(class_type,arg_type) \
	inline class_type::complex(const arg_type &x):ancestor(x) {}
	#define MATH_OPERATOR_DECL(class_type,arg_type,op) \
	/** \brief Operator op against arg_type. */ \
	class_type &operator op(const arg_type &);
	#define MATH_OPERATOR(class_type,arg_type,op) \
	inline class_type & class_type::operator op(const arg_type &other) \
	{ \
		return ancestor::operator op(other); \
	}
	#define COMPARISON_OPERATOR_DECL(class_type,arg_type) \
	/** \brief Operator op against arg_type. */ \
	bool operator==(const arg_type &) const;
	#define COMPARISON_OPERATOR(class_type,arg_type) \
	inline bool class_type::operator==(const arg_type &other) const \
	{ \
		return ancestor::operator==(other); \
	}

	// Forward declarations.
	template <>
	class complex<piranha::mp_rational>;
	template <>
	class complex<piranha::mp_integer>;

	/// Complex counterpart of piranha::mp_rational.
	/**
	 * Stores two piranha::mp_rational members as real and imaginary part.
	 */
	template <>
	class complex<piranha::mp_rational>:
		public piranha::complex_generic_mp_container<piranha::mp_rational,complex<piranha::mp_rational> >,
		boost::field_operators<complex<piranha::mp_rational>,
		boost::field_operators<complex<piranha::mp_rational>, int,
		boost::field_operators<complex<piranha::mp_rational>, double,
		boost::field_operators<complex<piranha::mp_rational>, piranha::mp_rational,
		boost::field_operators<complex<piranha::mp_rational>, piranha::mp_integer,
		boost::field_operators<complex<piranha::mp_rational>, complex<int>,
		boost::field_operators<complex<piranha::mp_rational>, complex<double>,
		boost::field_operators<complex<piranha::mp_rational>, complex<piranha::mp_integer>
		> > > > > > > >
	{
			typedef piranha::complex_generic_mp_container<piranha::mp_rational,complex<piranha::mp_rational> > ancestor;
		public:
			/// STL-like typedef for internal scalar type.
			typedef piranha::mp_rational value_type;
			/// Default constructor.
			explicit complex(): ancestor() {}
			/// Constructor from integer.
			explicit complex(const int &n): ancestor(n) {}
			/// Constructor from integer real and imaginary parts.
			explicit complex(const int &r, const int &i): ancestor(r,i) {}
			/// Constructor from complex integer.
			explicit complex(const complex<int> &c): ancestor(c) {}
			/// Constructor from double.
			explicit complex(const double &x): ancestor(x) {}
			/// Constructor from double real and imaginary parts..
			explicit complex(const double &r, const double &i): ancestor(r,i) {}
			/// Constructor from complex double.
			explicit complex(const complex<double> &c): ancestor(c) {}
			/// Constructor from value_type.
			explicit complex(const value_type &x): ancestor(x) {}
			/// Constructor from real and imaginary parts of type value_type.
			explicit complex(const value_type &r, const value_type &i): ancestor(r,i) {}
			/// Constructor from std::string.
			/**
			 * @see complex_generic_mp_container::complex_generic_mp_container(const std::string &)
			 */
			explicit complex(const string &str): ancestor(str) {}
			/// Constructor from C string.
			/**
			 * @see complex(const std::string &).
			 */
			explicit complex(const char *str): ancestor(str) {}
			MATH_OPERATOR_DECL(complex,int,=)
			MATH_OPERATOR_DECL(complex,int,+=)
			MATH_OPERATOR_DECL(complex,int,-=)
			MATH_OPERATOR_DECL(complex,int,*=)
			MATH_OPERATOR_DECL(complex,int,/=)
			COMPARISON_OPERATOR_DECL(complex,int)
			MATH_OPERATOR_DECL(complex,double,=)
			MATH_OPERATOR_DECL(complex,double,+=)
			MATH_OPERATOR_DECL(complex,double,-=)
			MATH_OPERATOR_DECL(complex,double,*=)
			MATH_OPERATOR_DECL(complex,double,/=)
			COMPARISON_OPERATOR_DECL(complex,double)
			MATH_OPERATOR_DECL(complex,complex,+=)
			MATH_OPERATOR_DECL(complex,complex,-=)
			MATH_OPERATOR_DECL(complex,complex,*=)
			MATH_OPERATOR_DECL(complex,complex,/=)
			COMPARISON_OPERATOR_DECL(complex,complex)
			MATH_OPERATOR_DECL(complex,complex<int>,=)
			MATH_OPERATOR_DECL(complex,complex<int>,+=)
			MATH_OPERATOR_DECL(complex,complex<int>,-=)
			MATH_OPERATOR_DECL(complex,complex<int>,*=)
			MATH_OPERATOR_DECL(complex,complex<int>,/=)
			COMPARISON_OPERATOR_DECL(complex,complex<int>)
			MATH_OPERATOR_DECL(complex,complex<double>,=)
			MATH_OPERATOR_DECL(complex,complex<double>,+=)
			MATH_OPERATOR_DECL(complex,complex<double>,-=)
			MATH_OPERATOR_DECL(complex,complex<double>,*=)
			MATH_OPERATOR_DECL(complex,complex<double>,/=)
			COMPARISON_OPERATOR_DECL(complex,complex<double>)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,+=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,-=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,*=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,/=)
			COMPARISON_OPERATOR_DECL(complex,piranha::mp_rational)
			// Interop with other mp types.
			CTOR_DECL(complex,piranha::mp_integer)
			CTOR_DECL(complex,complex<piranha::mp_integer>)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,+=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,-=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,*=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,/=)
			COMPARISON_OPERATOR_DECL(complex,piranha::mp_integer)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_integer>,=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_integer>,+=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_integer>,-=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_integer>,*=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_integer>,/=)
			COMPARISON_OPERATOR_DECL(complex,complex<piranha::mp_integer>)
			/// Exponentiation.
			/**
			 * @throws zero_division_error if power is negative and value is zero.
			 * @throws value_error if power cannot be calculated exactly.
			 */
			complex pow(const double &y) const
			{
				if (piranha::utils::is_integer(y)) {
					return pow_int((int)y);
				} else {
				}	return pow_double(y);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n is zero.
			 * @throws value_error if root cannot be calculated exactly.
			 */
			complex root(const int &n) const
			{
				return pow(1. / double(n));
			}
		private:
			complex pow_int(const int &n) const
			{
				complex retval;
				// For negative powers, we must guard against division by zero.
				if (n < 0) {
					if ((*this) == 0) {
						piranha_throw(zero_division_error,"cannot raise zero to negative power");
					} else {
						// If source is non-zero, we can invert it and the calculate the power simply by multiplying.
						retval = invert();
						const size_t count = static_cast<size_t>(-n);
						complex tmp(retval);
						for (size_t i = 1; i < count; ++i) {
							retval *= tmp;
						}
					}
				} else {
					retval = 1;
					const size_t count = static_cast<size_t>(n);
					for (size_t i = 0; i < count; ++i) {
						retval *= (*this);
					}
				}
				return retval;
			}
			complex pow_double(const double &y) const
			{
				complex retval;
				if (y < 0) {
					if ((*this) == 0) {
						piranha_throw(zero_division_error,"cannot raise zero to negative power");
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						piranha_throw(value_error,"cannot raise complex rational "
							"different from unity to negative real power");
					}
				} else if (y == 0) {
					// If y == 0, then x**0 == 1 for every x.
					retval = 1;
				} else {
					// If y > 0, we can accept only 0^y and 1^y.
					if ((*this) == 0) {
						retval = 0;
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						piranha_throw(value_error,"cannot raise complex rational "
							"different from unity to positive real power");
					}
				}
				return retval;
			}
	};

	MATH_OPERATOR(complex<piranha::mp_rational>,int,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,int,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,int,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,int,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,int,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,int)
	MATH_OPERATOR(complex<piranha::mp_rational>,double,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,double,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,double,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,double,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,double,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,double)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_rational>,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_rational>,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_rational>,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_rational>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_rational>)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<int>,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<int>,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<int>,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<int>,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<int>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,complex<int>)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<double>,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<double>,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<double>,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<double>,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<double>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,complex<double>)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_rational,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_rational,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_rational,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_rational,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_rational,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,piranha::mp_rational)

	/// Overload of standard swap function for std::complex<piranha::mp_rational>.
	/**
	 * Will use the swap() method internally.
	 */
	inline void swap(complex<piranha::mp_rational> &qc1, complex<piranha::mp_rational> &qc2)
	{
		qc1.swap(qc2);
	}

	/// Overload standard power function for std::complex<piranha::mp_rational> and double argument.
	/**
	 * @see std::complex<piranha::mp_rational>::pow.
	 */
	inline complex<piranha::mp_rational> pow(const complex<piranha::mp_rational> &qc, const double &y)
	{
		return qc.pow(y);
	}

	/// Overload standard power function for std::complex<piranha::mp_rational> and int argument.
	/**
	 * @see std::complex<piranha::mp_rational>::pow.
	 */
	inline complex<piranha::mp_rational> pow(const complex<piranha::mp_rational> &qc, const int &y)
	{
		return qc.pow(y);
	}

	/// Overload in stream operator>> for std::complex<piranha::mp_rational>.
	inline std::istream &operator>>(std::istream &i, complex<piranha::mp_rational> &qc)
	{
		string tmp_str;
		getline(i,tmp_str);
		complex<piranha::mp_rational> tmp(tmp_str);
		swap(qc,tmp);
		return i;
	}

	/// Overload hash_value function for std::complex<piranha::mp_rational>.
	inline size_t hash_value(const complex<piranha::mp_rational> &qc)
	{
		return qc.hash();
	}

	/// Complex counterpart of piranha::mp_integer.
	/**
	 * Stores two piranha::mp_integer members as real and imaginary part.
	 */
	template <>
	class complex<piranha::mp_integer>:
		public piranha::complex_generic_mp_container<piranha::mp_integer,complex<piranha::mp_integer> >,
		boost::field_operators<complex<piranha::mp_integer>,
		boost::field_operators<complex<piranha::mp_integer>, int,
		boost::field_operators<complex<piranha::mp_integer>, double,
		boost::field_operators<complex<piranha::mp_integer>, piranha::mp_integer,
		boost::field_operators<complex<piranha::mp_integer>, complex<int>,
		boost::field_operators<complex<piranha::mp_integer>, complex<double>
		> > > > > >
	{
			typedef piranha::complex_generic_mp_container<piranha::mp_integer,complex<piranha::mp_integer> > ancestor;
		public:
			/// STL-like typedef for internal scalar type.
			typedef piranha::mp_integer value_type;
			/// Default constructor.
			explicit complex(): ancestor() {}
			/// Constructor from integer.
			explicit complex(const int &n): ancestor(n) {}
			/// Constructor from integer real and imaginary parts.
			explicit complex(const int &r, const int &i): ancestor(r,i) {}
			/// Constructor from complex integer.
			explicit complex(const complex<int> &c): ancestor(c) {}
			/// Constructor from double.
			explicit complex(const double &x): ancestor(x) {}
			/// Constructor from double real and imaginary parts..
			explicit complex(const double &r, const double &i): ancestor(r,i) {}
			/// Constructor from complex double.
			explicit complex(const complex<double> &c): ancestor(c) {}
			/// Constructor from value_type.
			explicit complex(const value_type &x): ancestor(x) {}
			/// Constructor from real and imaginary parts of type value_type.
			explicit complex(const value_type &r, const value_type &i): ancestor(r,i) {}
			/// Constructor from std::string.
			/**
			 * @see complex_generic_mp_container::complex_generic_mp_container(const std::string &)
			 */
			explicit complex(const string &str): ancestor(str) {}
			/// Constructor from C string.
			/**
			 * @see complex(const std::string &).
			 */
			explicit complex(const char *str): ancestor(str) {}
			MATH_OPERATOR_DECL(complex,int,=)
			MATH_OPERATOR_DECL(complex,int,+=)
			MATH_OPERATOR_DECL(complex,int,-=)
			MATH_OPERATOR_DECL(complex,int,*=)
			MATH_OPERATOR_DECL(complex,int,/=)
			COMPARISON_OPERATOR_DECL(complex,int)
			MATH_OPERATOR_DECL(complex,double,=)
			MATH_OPERATOR_DECL(complex,double,+=)
			MATH_OPERATOR_DECL(complex,double,-=)
			MATH_OPERATOR_DECL(complex,double,*=)
			MATH_OPERATOR_DECL(complex,double,/=)
			COMPARISON_OPERATOR_DECL(complex,double)
			MATH_OPERATOR_DECL(complex,complex,+=)
			MATH_OPERATOR_DECL(complex,complex,-=)
			MATH_OPERATOR_DECL(complex,complex,*=)
			MATH_OPERATOR_DECL(complex,complex,/=)
			COMPARISON_OPERATOR_DECL(complex,complex)
			MATH_OPERATOR_DECL(complex,complex<int>,=)
			MATH_OPERATOR_DECL(complex,complex<int>,+=)
			MATH_OPERATOR_DECL(complex,complex<int>,-=)
			MATH_OPERATOR_DECL(complex,complex<int>,*=)
			MATH_OPERATOR_DECL(complex,complex<int>,/=)
			COMPARISON_OPERATOR_DECL(complex,complex<int>)
			MATH_OPERATOR_DECL(complex,complex<double>,=)
			MATH_OPERATOR_DECL(complex,complex<double>,+=)
			MATH_OPERATOR_DECL(complex,complex<double>,-=)
			MATH_OPERATOR_DECL(complex,complex<double>,*=)
			MATH_OPERATOR_DECL(complex,complex<double>,/=)
			COMPARISON_OPERATOR_DECL(complex,complex<double>)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,+=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,-=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,*=)
			MATH_OPERATOR_DECL(complex,piranha::mp_integer,/=)
			COMPARISON_OPERATOR_DECL(complex,piranha::mp_integer)
			// Interop with other mp types.
			CTOR_DECL(complex,piranha::mp_rational)
			CTOR_DECL(complex,complex<piranha::mp_rational>)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,+=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,-=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,*=)
			MATH_OPERATOR_DECL(complex,piranha::mp_rational,/=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_rational>,=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_rational>,+=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_rational>,-=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_rational>,*=)
			MATH_OPERATOR_DECL(complex,complex<piranha::mp_rational>,/=)
			/// Exponentiation.
			/**
			 * @throws zero_division_error if power is negative and value is zero.
			 * @throws value_error if power cannot be calculated exactly.
			 */
			complex pow(const double &y) const
			{
				if (piranha::utils::is_integer(y)) {
					return pow_int((int)y);
				} else {
				}	return pow_double(y);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n is zero.
			 * @throws value_error if root cannot be calculated exactly.
			 */
			complex root(const int &n) const
			{
				return pow(1. / double(n));
			}
		private:
			complex pow_int(const int &n) const
			{
				complex retval;
				// For negative powers, we must guard against division by zero.
				if (n < 0) {
					if ((*this) == 0) {
						piranha_throw(zero_division_error,"cannot divide by zero");
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						piranha_throw(value_error,"cannot raise complex integer number different from unity to "
							"negative integer power");
					}
				} else {
					retval = 1;
					const size_t count = static_cast<size_t>(n);
					for (size_t i = 0; i < count; ++i) {
						retval *= (*this);
					}
				}
				return retval;
			}
			complex pow_double(const double &y) const
			{
				complex retval;
				if (y < 0) {
					if ((*this) == 0) {
						piranha_throw(zero_division_error,"cannot divide by zero");
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						piranha_throw(value_error,"cannot raise complex integer number different from unity to "
							"negative real power");
					}
				} else if (y == 0) {
					// If y == 0, then x**0 == 1 for every x.
					retval = 1;
				} else {
					// If y > 0, we can accept only 0^y and 1^y.
					if ((*this) == 0) {
						retval = 0;
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						piranha_throw(value_error,"cannot raise complex integer different from unity to "
							"positive real power");
					}
				}
				return retval;
			}
	};

	MATH_OPERATOR(complex<piranha::mp_integer>,int,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,int,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,int,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,int,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,int,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,int)
	MATH_OPERATOR(complex<piranha::mp_integer>,double,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,double,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,double,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,double,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,double,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,double)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_integer>,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_integer>,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_integer>,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_integer>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_integer>)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<int>,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<int>,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<int>,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<int>,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<int>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,complex<int>)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<double>,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<double>,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<double>,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<double>,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<double>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,complex<double>)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_integer,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_integer,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_integer,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_integer,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_integer,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,piranha::mp_integer)

	// Mixed type operators.
	// Complex rational.
	CTOR(complex<piranha::mp_rational>,piranha::mp_integer)
	CTOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_integer,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_integer,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_integer,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_integer,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,piranha::mp_integer,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,piranha::mp_integer)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>,=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>,+=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>,-=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>,*=)
	MATH_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>,/=)
	COMPARISON_OPERATOR(complex<piranha::mp_rational>,complex<piranha::mp_integer>)
	// Complex integer.
	CTOR(complex<piranha::mp_integer>,piranha::mp_rational)
	CTOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_rational,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_rational,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_rational,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_rational,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,piranha::mp_rational,/=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,/=)

	/// Overload of standard swap function for std::complex<piranha::mp_integer>.
	/**
	 * Will use the swap() method internally.
	 */
	inline void swap(complex<piranha::mp_integer> &zc1, complex<piranha::mp_integer> &zc2)
	{
		zc1.swap(zc2);
	}

	/// Overload standard power function for std::complex<piranha::mp_integer> and double argument.
	/**
	 * @see std::complex<piranha::mp_integer>::pow.
	 */
	inline complex<piranha::mp_integer> pow(const complex<piranha::mp_integer> &zc, const double &y)
	{
		return zc.pow(y);
	}

	/// Overload standard power function for std::complex<piranha::mp_integer> and int argument.
	/**
	 * @see std::complex<piranha::mp_integer>::pow.
	 */
	inline complex<piranha::mp_integer> pow(const complex<piranha::mp_integer> &zc, const int &y)
	{
		return zc.pow(y);
	}

	/// Overload in stream operator>> for std::complex<piranha::mp_integer>.
	inline std::istream &operator>>(std::istream &i, complex<piranha::mp_integer> &zc)
	{
		string tmp_str;
		getline(i,tmp_str);
		complex<piranha::mp_integer> tmp(tmp_str);
		swap(zc,tmp);
		return i;
	}

	/// Overload hash_value function for std::complex<piranha::mp_integer>.
	inline size_t hash_value(const complex<piranha::mp_integer> &zc)
	{
		return zc.hash();
	}

	#undef CTOR_DECL
	#undef CTOR
	#undef MATH_OPERATOR_DECL
	#undef MATH_OPERATOR
	#undef COMPARISON_OPERATOR_DECL
	#undef COMPARISON_OPERATOR
}

#endif
