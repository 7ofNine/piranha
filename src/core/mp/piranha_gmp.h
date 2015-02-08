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

#ifndef PIRANHA_PIRANHA_GMP_H
#define PIRANHA_PIRANHA_GMP_H

#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/operators.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../math.h"
#include "complex_generic_mp_container.h"
#include "mp_commons.h"

// TODO:
// - fix wrong behaviour in Python whe constructing from big ints that get converted to double... this is to be done from the Python exposition by overloading __init__.
// - better performance for complex ints using multadd (possibly through stl-like functor) and completion of API:
// - better handling when building mp_integer from bogus string
// - full interaction with POD types. Should not be much too effort now to complete...
// - do not use boost::operator ++ and --, use common impl in mp_commons instead and overload correctly for pre/post increment.
// - overload std::pow for exponentiation to rational
// - interaction between real of one type with complex of other type?
// - better pow for complex, like handling the case in which there is only real or imaginary part and pow can be forwarded
//   to the scalar implementation? Maybe overkill...
// ...

namespace piranha
{
	// Forward declaration of classes.
	class mp_rational;
	class mp_integer;
}

namespace std
{
	// Forward declarations.
	template <>
	class complex<piranha::mp_rational>;
	template <>
	class complex<piranha::mp_integer>;
}

namespace piranha
{
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
		normal_check(other); \
		m_value op other __VA_ARGS__; \
		return *this; \
	}
	#define FORWARDING_DIVISION_OPERATOR_DECL(class_type,arg_type) \
	/** \brief Operator /= against arg_type. */ \
	class_type & operator /=(const arg_type &);
	#define FORWARDING_DIVISION_OPERATOR(class_type,arg_type,...) \
	inline class_type & class_type::operator /=(const arg_type &other) \
	{ \
		normal_check(other); \
		if (other == 0) { \
			PIRANHA_THROW(zero_division_error,"cannot divide by zero"); \
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
	#define COMPLEX_COMPARISON_OPERATOR(arg_type) \
	/** \brief Comparison against std::complex< arg_type >. */ \
	bool operator==(const std::complex< arg_type > &c) const\
	{ \
		return complex_comparison(c);\
	}

	/// Multiprecision rational class.
	/**
	 * Wraps a GMP mpq_class.
	 */
	class mp_rational:
		public gmp_toolbox<mpq_class,mp_rational>,
		boost::incrementable<mp_rational,
		boost::decrementable<mp_rational,
		boost::ordered_field_operators<mp_rational,
		boost::ordered_field_operators<mp_rational, int,
		boost::ordered_field_operators<mp_rational, double,
		boost::ordered_field_operators<mp_rational, mp_integer,
		boost::equality_comparable<mp_rational, std::complex<double>,
		boost::equality_comparable<mp_rational, std::complex<int>
		> > > > > > > >
	{
			friend class gmp_toolbox<mpq_class,mp_rational>;
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
			/// Constructor from int.
			explicit mp_rational(const int &n):m_value(n) {}
			/// Constructor from integer numerator and denominator.
			/**
			 * @throws zero_division_error if denominator is zero.
			 */
			explicit mp_rational(const int &n, const int &d):m_value(0)
			{
				construct_from_numden(n,d);
			}
			/// Constructor from double.
			explicit mp_rational(const double &x):m_value(0)
			{
				normal_check(x);
				m_value = x;
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
			COMPLEX_COMPARISON_OPERATOR(double)
			COMPLEX_COMPARISON_OPERATOR(int)
			// Interoperability with mp_integer.
			/// Constructor from piranha::mp_integer numerator and denominator.
			/**
			 * @throws zero_division_error if denominator is zero.
			 */
			explicit mp_rational(const mp_integer &, const mp_integer &);
			FORWARDING_CTOR_DECL(mp_rational,mp_integer)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_rational,mp_integer,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_rational,mp_integer)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_integer,==)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_integer,<)
			FORWARDING_COMPARISON_OPERATOR_DECL(mp_rational,mp_integer,>)
			/// Get copy of numerator.
			mp_integer get_num() const;
			/// Get copy of denominator.
			mp_integer get_den() const;
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
			/// Convert to int.
			/**
			 * @throws value_error if the denominator is not unitary.
			 * @throws std::overflow_error if the numerator overflows int type.
			 */
			int to_int() const
			{
				if (m_value.get_den() != 1) {
					PIRANHA_THROW(value_error,"cannot convert rational to int if denominator is non-unitary");
				}
				// NOTE: here we test for fits into signed integer, but the return value below will be long signed
				// integer because that's what offered by the API. It will be safely downcast to int on exit.
				if (!m_value.get_num().fits_sint_p()) {
					PIRANHA_THROW(std::overflow_error,"numerator is too large while converting rational to int");
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
			std::size_t hash() const
			{
				std::size_t retval = 0;
				const __mpz_struct *num = mpq_numref(m_value.get_mpq_t()), *den = mpq_denref(m_value.get_mpq_t());
				const std::size_t num_limb_size = std::abs(num->_mp_size), den_limb_size = std::abs(den->_mp_size);
				for (std::size_t i = 0; i < num_limb_size; ++i) {
					boost::hash_combine(retval, num->_mp_d[i]);
				}
				for (std::size_t i = 0; i < den_limb_size; ++i) {
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
			/// Binomial coefficient of this over k.
			/**
			 * Internally it will use the piranha::generic_choose function.
			 * @param[in] k integer over which binomial coefficient of this is calculated.
			 * @param[out] retval binomial coefficient of this over k.
			 */
			mp_rational choose(const int &k) const
			{
				return generic_choose(*this,k);
			}
			/// Exponentiation.
			/**
			 * @throws zero_division_error if power is negative and value is zero.
			 * @throws value_error if power cannot be calculated exactly.
			 */
			mp_rational pow(const double &y) const
			{
				normal_check(y);
				if (is_integer(y)) {
					return pow_int((int)y);
				} else {
					return pow_double(y);
				}
			}
			/// Rational exponentiation.
			mp_rational pow(const mp_rational &q) const
			{
				return rat_pow(*this,q);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n_ is zero.
			 * @throws value_error if root cannot be calculated exactly or this is negative.
			 */
			mp_rational root(const int &n_) const
			{
				if (m_value < 0 && n_ != 1 && n_ != -1) {
					PIRANHA_THROW(value_error,"cannot calculate root of negative rational number");
				}
				mp_rational retval;
				if (n_ == 0) {
					PIRANHA_THROW(zero_division_error,"cannot calculate zero-th root of rational number");
				} else if (n_ == 1) {
					retval = *this;
					return retval;
				} else if (n_ == -1) {
					retval = pow(-1);
					return retval;
				}
				const std::size_t n = (n_ > 0) ? n_ : -n_;
				if (!mpz_root(mpq_numref(retval.m_value.get_mpq_t()),mpq_numref(m_value.get_mpq_t()),n) ||
					!mpz_root(mpq_denref(retval.m_value.get_mpq_t()),mpq_denref(m_value.get_mpq_t()),n)) {
					PIRANHA_THROW(value_error,"rational number is not an exact nth root");
				}
				// Better to canonicalise, for peace of mind.
				retval.m_value.canonicalize();
				if (n_ < 0) {
					// Let's guard against division by zero below.
					if (retval == 0) {
						PIRANHA_THROW(zero_division_error,"cannot calculate negative root of zero");
					}
					mpq_inv(retval.m_value.get_mpq_t(), mpq_class(retval.m_value).get_mpq_t());
				}
				return retval;
			}
			/// Print to stream in Tex format.
			void printTex(std::ostream &outStream) const
			{
				if (m_value.get_den() == 1) {
					outStream << ' ' << m_value.get_num() << ' ';
				} else {
					if (m_value.get_num() >= 0) {
						outStream << "\\frac{" << m_value.get_num() << "}{" << m_value.get_den() << '}';
					} else {
						outStream << "-\\frac{" << (-m_value.get_num()) << "}{" << m_value.get_den() << '}';
					}
				}
			}
		private:
			template <class T>
			void construct_from_numden(const T &n, const T &d)
			{
				// Guard against division by zero.
				if (d == 0) {
					PIRANHA_THROW(zero_division_error,"cannot create rational with zero as denominator");
				}
				mpq_class tmp(n,d);
				mpq_canonicalize(tmp.get_mpq_t());
				// Swap content (more efficient than copying).
				mpz_swap(mpq_numref(m_value.get_mpq_t()),mpq_numref(tmp.get_mpq_t()));
				mpz_swap(mpq_denref(m_value.get_mpq_t()),mpq_denref(tmp.get_mpq_t()));
			}
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
							PIRANHA_THROW(value_error,"invalid string input");
						}
						break;
					case 2:
						{
						mpz_class num(0), den(0);
						try {
							num = mpz_class(split_v[0]);
							den = mpz_class(split_v[1]);
						} catch (const std::invalid_argument &) {
							PIRANHA_THROW(value_error,"invalid string input");
						}
						if (den == 0) {
							PIRANHA_THROW(zero_division_error,"cannot create rational with zero as denominator");
						}
						mpq_class tmp_q(num,den);
						mpq_canonicalize(tmp_q.get_mpq_t());
						m_value = tmp_q;
						}
						break;
					default:
						PIRANHA_THROW(value_error,"invalid string input");
				}
			}
			mp_rational pow_int(const int &n) const {
				mp_rational retval;
				if (m_value == 0 && n < 0) {
					PIRANHA_THROW(zero_division_error,"cannot raise zero to negative power");
				}
				if (n < 0) {
					mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (unsigned long)(-n));
					mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (unsigned long)(-n));
					// We need to canonicalize, since negative numbers may have gone to the denominator.
					retval.m_value.canonicalize();
				} else {
					mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (unsigned long)n);
					mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (unsigned long)n);
				}
				return retval;
			}
			mp_rational pow_double(const double &y) const {
				mp_rational retval;
				// If negative, only 1^-something is reasonable.
				if (y < 0) {
					if (m_value == 0) {
						PIRANHA_THROW(zero_division_error,"cannot raise zero to negative power");
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						PIRANHA_THROW(value_error,"cannot raise rational number different from unity to "
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
						PIRANHA_THROW(value_error,"cannot raise rational number different from unity to "
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

	/// Multiprecision integer class.
	/**
	 * Wraps a GMP mpz_class.
	 */
	class mp_integer:
		public gmp_toolbox<mpz_class,mp_integer>,
		boost::incrementable<mp_integer,
		boost::decrementable<mp_integer,
		boost::ordered_field_operators<mp_integer,
		boost::ordered_field_operators<mp_integer, int,
		boost::ordered_field_operators<mp_integer, double,
		boost::modable<mp_integer,
		boost::modable<mp_integer, int,
		boost::equality_comparable<mp_integer, std::complex<double>,
		boost::equality_comparable<mp_integer, std::complex<int>
		> > > > > > > > >
	{
			friend class gmp_toolbox<mpz_class,mp_integer>;
		public:
			/// Make friends with piranha::mp_rational.
			friend class mp_rational;
			/// Default constructor.
			/**
			 * Initialises value to zero.
			 */
			explicit mp_integer():m_value(0) {}
			/// Constructor from std::string.
			explicit mp_integer(const std::string &str):m_value(0)
			{
				try {
					m_value = mpz_class(str);
				} catch (const std::invalid_argument &) {
					PIRANHA_THROW(value_error,"invalid string input");
				}
			}
			/// Constructor from C string.
			explicit mp_integer(const char *str):m_value(str) {}
			/// Constructor from int.
			explicit mp_integer(const int &n):m_value(n) {}
			/// Constructor from double.
			explicit mp_integer(const double &x):m_value(0)
			{
				normal_check(x);
				m_value = x;
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
			COMPLEX_COMPARISON_OPERATOR(double)
			COMPLEX_COMPARISON_OPERATOR(int)
			// Interoperability with mp_rational.
			FORWARDING_CTOR_DECL(mp_integer,mp_rational)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,+=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,-=)
			FORWARDING_MATH_OPERATOR_DECL(mp_integer,mp_rational,*=)
			FORWARDING_DIVISION_OPERATOR_DECL(mp_integer,mp_rational)
			/// Mod operator.
			mp_integer &operator%=(const int &n)
			{
				m_value %= n;
				return *this;
			}
			/// Mod operator.
			mp_integer &operator%=(const mp_integer &z)
			{
				m_value %= z.m_value;
				return *this;
			}
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
			/// Convert to int.
			/**
			 * @throws std::overflow_error if the numerator overflows int type.
			 */
			int to_int() const
			{
				if (!m_value.fits_sint_p()) {
					PIRANHA_THROW(std::overflow_error,"multiprecision integer too big to be converted to int");
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
			/// Multiply-accumulate.
			/**
			 * Will use GMP's addmul() function, see http://gmplib.org/manual/Integer-Arithmetic.html#Integer-Arithmetic.
			 */
			void multiply_accumulate(const mp_integer &y, const mp_integer &z)
			{
				mpz_addmul(m_value.get_mpz_t(),y.m_value.get_mpz_t(),z.m_value.get_mpz_t());
			}
			/// Least common multiplier.
			/**
			 * Set this to the least common multiplier of a and b.
			 * Will use GMP's mpz_lcm() function, see
			 * http://gmplib.org/manual/Number-Theoretic-Functions.html#Number-Theoretic-Functions.
			 */
			void lcm(const mp_integer &a, const mp_integer &b)
			{
				mpz_lcm(m_value.get_mpz_t(),a.m_value.get_mpz_t(),b.m_value.get_mpz_t());
			}
			/// Hash value.
			/**
			 * Internally uses boost::hash_combine on the GMP limbs of the number.
			 */
			std::size_t hash() const
			{
				std::size_t retval = 0;
				const __mpz_struct *ptr = m_value.get_mpz_t();
				const std::size_t limb_size = std::abs(ptr->_mp_size);
				for (std::size_t i = 0; i < limb_size; ++i) {
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
			/// Factorial.
			/**
			 * Internally it will use the GMP mpz_fac_ui function.
			 * @param[out] retval factorial of this.
			 * @throws value_error if this is negative.
			 * @throws std::overflow_error if this overflows int type.
			 */
			mp_integer factorial() const
			{
				if ((*this) < 0) {
					PIRANHA_THROW(value_error,"cannot calculate factorial of negative integer");
				}
				mp_integer retval;
				mpz_fac_ui(retval.m_value.get_mpz_t(),(unsigned long)to_int());
				return retval;
			}
			/// Double factorial.
			/**
			 * Internally it will use the piranha::generic_double_factorial function.
			 * @param[out] retval double factorial of this.
			 * @see piranha::generic_double_factorial.
			 */
			mp_integer double_factorial() const
			{
				return generic_double_factorial(*this);
			}
			/// Binomial coefficient of this over k.
			/**
			 * Internally it will use the GMP mpz_bin_ui function.
			 * @param[in] k integer over which binomial coefficient of this is calculated.
			 * @param[out] retval binomial coefficient of this over k.
			 */
			mp_integer choose(const int &k) const
			{
				mp_integer retval(0);
				if (k < 0) {
					return retval;
				}
				mpz_bin_ui(retval.m_value.get_mpz_t(),m_value.get_mpz_t(),(unsigned long)k);
				return retval;
			}
			/// Exponentiation.
			/**
			 * @throws zero_division_error if power is negative and value is zero.
			 * @throws value_error if power cannot be calculated exactly.
			 */
			mp_integer pow(const double &y) const
			{
				normal_check(y);
				if (is_integer(y)) {
					return pow_int((int)y);
				} else {
					return pow_double(y);
				}
			}
			/// Rational exponentiation.
			mp_integer pow(const mp_rational &q) const
			{
				return rat_pow(*this,q);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n_ is zero.
			 * @throws value_error if root cannot be calculated exactly or this is negative.
			 */
			mp_integer root(const int &n_) const
			{
				if (m_value < 0) {
					PIRANHA_THROW(value_error,"cannot calculate nth-root of negative integer number");
				}
				mp_integer retval;
				if (n_ == 0) {
					PIRANHA_THROW(zero_division_error,"cannot calculate zero-th root");
				} else if (n_ == 1) {
					retval = *this;
					return retval;
				} else if (n_ < 0) {
					PIRANHA_THROW(value_error,"integer numbers different from unity cannot be arguments of negative root");
				}
				const std::size_t n = static_cast<std::size_t>(n_);
				if (!mpz_root(retval.m_value.get_mpz_t(),m_value.get_mpz_t(),n)) {
					PIRANHA_THROW(value_error,"integer coefficient is not an exact nth root");
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
						PIRANHA_THROW(zero_division_error,"cannot divide by zero");
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						PIRANHA_THROW(value_error,"cannot raise integer number different from unity to "
							"negative integer power");
					}
				} else {
					mpz_pow_ui(retval.m_value.get_mpz_t(), m_value.get_mpz_t(), (std::size_t)n);
				}
				return retval;
			}
			mp_integer pow_double(const double &y) const
			{
				mp_integer retval;
				// If negative, only 1^-something is reasonable.
				if (y < 0) {
					if (m_value == 0) {
						PIRANHA_THROW(zero_division_error,"cannot divide by zero");
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						PIRANHA_THROW(value_error,"cannot raise integer number different from unity to negative real power");
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
						PIRANHA_THROW(value_error,"cannot raise integer number different from unity to positive real power");
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

	// Mixed operations between mp types.
	inline mp_rational::mp_rational(const mp_integer &n, const mp_integer &d):m_value(0)
	{
		construct_from_numden(n.get_internal(),d.get_internal());
	}
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

	// mp_rational's operations that require the definition of mp_integer.
	inline mp_integer mp_rational::get_num() const
	{
		mp_integer retval;
		retval.m_value = m_value.get_num();
		return retval;
	}

	inline mp_integer mp_rational::get_den() const
	{
		mp_integer retval;
		retval.m_value = m_value.get_den();
		return retval;
	}

	#undef FORWARDING_CTOR_DECL
	#undef FORWARDING_CTOR
	#undef FORWARDING_MATH_OPERATOR_DECL
	#undef FORWARDING_MATH_OPERATOR
	#undef FORWARDING_DIVISION_OPERATOR_DECL
	#undef FORWARDING_DIVISION_OPERATOR
	#undef FORWARDING_COMPARISON_OPERATOR_DECL
	#undef FORWARDING_COMPARISON_OPERATOR
	#undef COMPLEX_COMPARISON_OPERATOR
}

namespace std
{
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
				piranha::normal_check(y);
				if (piranha::is_integer(y)) {
					return pow_int((int)y);
				} else {
					return pow_double(y);
				}
			}
			/// Rational exponentiation.
			complex pow(const piranha::mp_rational &q) const
			{
				return rat_pow(*this,q);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n is zero.
			 * @throws value_error if root cannot be calculated exactly.
			 */
			complex root(const int &n) const
			{
				if (!n) {
					PIRANHA_THROW(zero_division_error,"cannot calculate zero-th root of complex rational number");
				}
				return pow(1. / double(n));
			}
		private:
			complex pow_int(const int &n) const
			{
				complex retval;
				// For negative powers, we must guard against division by zero.
				if (n < 0) {
					if ((*this) == 0) {
						PIRANHA_THROW(zero_division_error,"cannot raise zero to negative power");
					} else {
						// If source is non-zero, we can invert it and the calculate the power simply by multiplying.
						retval = invert();
						const std::size_t count = static_cast<std::size_t>(-n);
						complex tmp(retval);
						for (std::size_t i = 1; i < count; ++i) {
							retval *= tmp;
						}
					}
				} else {
					retval = 1;
					const std::size_t count = static_cast<std::size_t>(n);
					for (std::size_t i = 0; i < count; ++i) {
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
						PIRANHA_THROW(zero_division_error,"cannot raise zero to negative power");
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						PIRANHA_THROW(value_error,"cannot raise complex rational "
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
						PIRANHA_THROW(value_error,"cannot raise complex rational "
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
		boost::field_operators<complex<piranha::mp_integer>, piranha::mp_rational,
		boost::field_operators<complex<piranha::mp_integer>, complex<int>,
		boost::field_operators<complex<piranha::mp_integer>, complex<double>
		> > > > > > >
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
			COMPARISON_OPERATOR_DECL(complex,piranha::mp_rational)
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
				piranha::normal_check(y);
				if (piranha::is_integer(y)) {
					return pow_int((int)y);
				} else {
					return pow_double(y);
				}
			}
			/// Rational exponentiation.
			complex pow(const piranha::mp_rational &q) const
			{
				return rat_pow(*this,q);
			}
			/// N-th root.
			/**
			 * @throws zero_division_error if n is zero.
			 * @throws value_error if root cannot be calculated exactly.
			 */
			complex root(const int &n) const
			{
				if (!n) {
					PIRANHA_THROW(zero_division_error,"cannot calculate zero-th root of complex integer number");
				}
				return pow(1. / double(n));
			}
		private:
			complex pow_int(const int &n) const
			{
				complex retval;
				// For negative powers, we must guard against division by zero.
				if (n < 0) {
					if ((*this) == 0) {
						PIRANHA_THROW(zero_division_error,"cannot divide by zero");
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						PIRANHA_THROW(value_error,"cannot raise complex integer number different from unity to "
							"negative integer power");
					}
				} else {
					retval = 1;
					const std::size_t count = static_cast<std::size_t>(n);
					for (std::size_t i = 0; i < count; ++i) {
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
						PIRANHA_THROW(zero_division_error,"cannot divide by zero");
					} else if ((*this) == 1) {
						retval = 1;
					} else {
						PIRANHA_THROW(value_error,"cannot raise complex integer number different from unity to "
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
						PIRANHA_THROW(value_error,"cannot raise complex integer different from unity to "
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
	COMPARISON_OPERATOR(complex<piranha::mp_integer>,piranha::mp_rational)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,+=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,-=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,*=)
	MATH_OPERATOR(complex<piranha::mp_integer>,complex<piranha::mp_rational>,/=)

	#undef CTOR_DECL
	#undef CTOR
	#undef MATH_OPERATOR_DECL
	#undef MATH_OPERATOR
	#undef COMPARISON_OPERATOR_DECL
	#undef COMPARISON_OPERATOR
}

#endif
