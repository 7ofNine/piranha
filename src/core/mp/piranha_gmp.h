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

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/operators.hpp>
#include <cmath>
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"

namespace piranha
{
	// Construct hash value from mpq_class.
	static inline size_t hash_from_mpq_class(const mpq_class &q)
	{
		size_t retval = 0;
		const __mpz_struct *num = mpq_numref(q.get_mpq_t()), *den = mpq_denref(q.get_mpq_t());
		const size_t num_limb_size = std::abs(num->_mp_size), den_limb_size = std::abs(den->_mp_size);
		for (size_t i = 0; i < num_limb_size; ++i) {
			boost::hash_combine(retval, num->_mp_d[i]);
		}
		for (size_t i = 0; i < den_limb_size; ++i) {
			boost::hash_combine(retval, den->_mp_d[i]);
		}
		return retval;
	}

	// Macros for operators.
	#define EQUALITY_OPERATOR(type) \
	/** \brief Test equality with type. */ \
	bool operator==(const type &other) const \
	{ \
		return (m_value == other); \
	}
	#define COMPARISON_OPERATOR(type) \
	/** \brief Comparison operator with type. */ \
	bool operator<(const type &other) const \
	{ \
		return (m_value < other); \
	}
	#define INPLACE_OPERATOR(self,op,type) \
	/** \brief In-place operator op for type. */ \
	self &operator op(const type &other) \
	{ \
		m_value op other; \
		return *this; \
	}
	#define INPLACE_DIVISION(self,type) \
	/** \brief In-place division with type. */ \
	self &operator/=(const type &other) \
	{ \
		if (other == 0) { \
			piranha_throw(zero_division_error,"cannot divide by zero"); \
		} \
		m_value /= other; \
		return *this; \
	}

	/// Multiprecision rational class.
	/**
	 * Wraps a GMP mpq_class. Interoperability with C++ int and double types is provided through constructors and
	 * in-place mathematical operators. Non in-place operators are provided through inheritance from the Boost operator library.
	 */
	class mp_rational:
		boost::ordered_field_operators<mp_rational,
		boost::ordered_field_operators<mp_rational, int,
		boost::ordered_field_operators<mp_rational, double
		> > >
	{
			friend std::ostream &operator<<(std::ostream &, const mp_rational &);
			typedef mp_rational self;
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
			explicit mp_rational(const std::string &str):m_value()
			{
				construct_from_string(str.c_str());
			}
			/// Constructor from C string.
			/**
			 * @see mp_rational(const std::string &).
			 */
			explicit mp_rational(const char *str):m_value()
			{
				construct_from_string(str);
			}
			/// Constructor from integer.
			explicit mp_rational(const int &n):m_value(n) {}
			/// Constructor from integer numerator and denominator.
			/**
			 * @throws zero_division_error if denominator is zero.
			 */
			explicit mp_rational(const int &n, const int &d):m_value()
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
			/// Swap content.
			/**
			 * Internally uses the mpz_swap GMP function on numerator and denominator..
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
				return hash_from_mpq_class(m_value);
			}
			/// Equality operator.
			bool operator==(const mp_rational &other) const
			{
				return (m_value == other.m_value);
			}
			/// Comparison operator.
			bool operator<(const mp_rational &other) const
			{
				return (m_value < other.m_value);
			}
			/// In-place addition.
			mp_rational &operator+=(const mp_rational &other)
			{
				m_value += other.m_value;
				return *this;
			}
			/// In-place subtraction.
			mp_rational &operator-=(const mp_rational &other)
			{
				m_value -= other.m_value;
				return *this;
			}
			/// In-place multiplication.
			mp_rational &operator*=(const mp_rational &other)
			{
				m_value *= other.m_value;
				return *this;
			}
			/// In-place division.
			/**
			 * @throws zero_division_error if dividing by zero.
			 */
			mp_rational &operator/=(const mp_rational &other)
			{
				if (other.m_value == 0) {
					piranha_throw(zero_division_error,"cannot divide by zero");
				}
				m_value /= other.m_value;
				return *this;
			}
			/// Negate in-plcae.
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
			// Operators against types directly compatible with mpq_class.
			EQUALITY_OPERATOR(int)
			EQUALITY_OPERATOR(double)
			COMPARISON_OPERATOR(int)
			COMPARISON_OPERATOR(double)
			INPLACE_OPERATOR(mp_rational,=,double)
			INPLACE_OPERATOR(mp_rational,=,int)
			INPLACE_OPERATOR(mp_rational,+=,int)
			INPLACE_OPERATOR(mp_rational,+=,double)
			INPLACE_OPERATOR(mp_rational,-=,int)
			INPLACE_OPERATOR(mp_rational,-=,double)
			INPLACE_OPERATOR(mp_rational,*=,int)
			INPLACE_OPERATOR(mp_rational,*=,double)
			INPLACE_DIVISION(mp_rational,int)
			INPLACE_DIVISION(mp_rational,double)
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
			mpq_class m_value;
	};

	#undef EQUALITY_OPERATOR
	#undef COMPARISON_OPERATOR
	#undef INPLACE_OPERATOR
	#undef INPLACE_DIVISION

	/// Overload out stream operator<< for piranha::mp_rational.
	inline std::ostream &operator<<(std::ostream &o, const mp_rational &q)
	{
		o << q.m_value;
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

	// Useful macros for operators.
	#define INPLACE_REAL_OPERATOR(op,type) \
	/** \brief In-place operator op for type. */ \
	complex &operator op(const type &x) \
	{ \
		m_real op x; \
		return *this; \
	}
	#define INPLACE_COMPLEX_OPERATOR(op,type) \
	/** \brief In-place operator op for complex type. */ \
	complex &operator op(const complex< type > &c) \
	{ \
		m_real op c.real(); \
		m_imag op c.imag(); \
		return *this; \
	}
	#define INPLACE_COMPLEX_MULT(type) \
	/** \brief In-place operator *= for complex type. */ \
	complex &operator*=(const complex< type > &c) \
	{ \
		return mult_by_complex(c); \
	}
	#define INPLACE_COMPLEX_DIV(type) \
	/** \brief In-place operator /= for complex type. */ \
	complex &operator/=(const complex< type > &c) \
	{ \
		return divide_by_complex(c); \
	}
	#define REAL_EQUALITY(type) \
	/** \brief Equality operator for type. */ \
	bool operator==(const type &x) \
	{ \
		return (m_real == x && m_imag == 0); \
	}
	#define COMPLEX_EQUALITY(type) \
	/** \brief Equality operator for complex type. */ \
	bool operator==(const complex< type > &c) \
	{ \
		return (m_real == c.real() && m_imag == c.imag()); \
	}

	/// Complex counterpart of piranha::mp_rational.
	template <>
	class complex<piranha::mp_rational>:
		boost::field_operators<complex<piranha::mp_rational>,
		boost::field_operators<complex<piranha::mp_rational>, int,
		boost::field_operators<complex<piranha::mp_rational>, double,
		boost::field_operators<complex<piranha::mp_rational>, piranha::mp_rational
		> > > >
	{
			friend ostream &operator<<(ostream &, const complex &);
		public:
			/// STL-like typedef for internal scalar type.
			typedef piranha::mp_rational value_type;
			/// Default constructor.
			explicit complex(): m_real(),m_imag() {}
			/// Constructor from integer.
			explicit complex(const int &n): m_real(n),m_imag() {}
			/// Constructor from integer real and imaginary parts..
			explicit complex(const int &r, const int &i): m_real(r),m_imag(i) {}
			/// Constructor from complex integer.
			explicit complex(const complex<int> &c): m_real(c.real()),m_imag(c.imag()) {}
			/// Constructor from double.
			explicit complex(const double &x): m_real(x),m_imag() {}
			/// Constructor from double real and imaginary parts..
			explicit complex(const double &r, const double &i): m_real(r),m_imag(i) {}
			/// Constructor from complex double.
			explicit complex(const complex<double> &c): m_real(c.real()),m_imag(c.imag()) {}
			/// Constructor from value_type.
			explicit complex(const value_type &x): m_real(x),m_imag() {}
			/// Constructor from real and imaginary parts of type value_type.
			explicit complex(const value_type &r, const value_type &i): m_real(r),m_imag(i) {}
			/// Constructor from std::string.
			/**
			 * Will raise a value_error exception if string is not valid. A valid string is of the form
			 * "(real,imag)", where "real" and "imag" are strings that can construct successfully a
			 * piranha::mp_rational.
			 * @see piranha::mp_rational::mp_rational(const std::string &).
			 * @throws value_error if string is invalid.
			 * @throws zero_division_error if string is valid but at least one denominator is zero.
			 */
			explicit complex(const string &str): m_real(),m_imag()
			{
				construct_from_string(str.c_str());
			}
			/// Constructor from C string.
			/**
			 * @see complex(const std::string &).
			 */
			explicit complex(const char *str): m_real(),m_imag()
			{
				construct_from_string(str);
			}
			/// Get const reference to the real part.
			const value_type &real() const
			{
				return m_real;
			}
			/// Get const reference to the imaginary part.
			const value_type &imag() const
			{
				return m_imag;
			}
			/// Swap content.
			void swap(complex &other)
			{
				m_real.swap(other.m_real);
				m_imag.swap(other.m_imag);
			}
			/// Hash value.
			/**
			 * Uses boost::hash_combine to combine the hashes of real and imaginary parts.
			 */
			size_t hash() const
			{
				boost::hash<value_type> hasher;
				size_t retval = hasher(m_real);
				boost::hash_combine(retval, hasher(m_imag));
				return retval;
			}
			/// Equality operator.
			bool operator==(const complex &other) const
			{
				return (m_real == other.m_real && m_imag == other.m_imag);
			}
			// Math operators.
			/// In-place addition.
			complex &operator+=(const complex &other)
			{
				m_real += other.m_real;
				m_imag += other.m_imag;
				return *this;
			}
			/// In-place subtraction.
			complex &operator-=(const complex &other)
			{
				m_real -= other.m_real;
				m_imag -= other.m_imag;
				return *this;
			}
			/// In-place multiplication.
			complex &operator*=(const complex &other)
			{
				return mult_by_complex(other);
			}
			/// In-place division.
			complex &operator/=(const complex &other)
			{
				return divide_by_complex(other);
			}
			// Maths for other types.
			REAL_EQUALITY(int)
			REAL_EQUALITY(double)
			REAL_EQUALITY(value_type)
			COMPLEX_EQUALITY(int)
			COMPLEX_EQUALITY(double)
			INPLACE_REAL_OPERATOR(=,int)
			INPLACE_REAL_OPERATOR(=,double)
			INPLACE_REAL_OPERATOR(=,value_type)
			INPLACE_COMPLEX_OPERATOR(=,int)
			INPLACE_COMPLEX_OPERATOR(=,double)
			INPLACE_REAL_OPERATOR(+=,int)
			INPLACE_REAL_OPERATOR(+=,double)
			INPLACE_REAL_OPERATOR(+=,value_type)
			INPLACE_COMPLEX_OPERATOR(+=,int)
			INPLACE_COMPLEX_OPERATOR(+=,double)
			INPLACE_REAL_OPERATOR(-=,int)
			INPLACE_REAL_OPERATOR(-=,double)
			INPLACE_REAL_OPERATOR(-=,value_type)
			INPLACE_COMPLEX_OPERATOR(-=,int)
			INPLACE_COMPLEX_OPERATOR(-=,double)
			INPLACE_REAL_OPERATOR(*=,int)
			INPLACE_REAL_OPERATOR(*=,double)
			INPLACE_REAL_OPERATOR(*=,value_type)
			INPLACE_REAL_OPERATOR(/=,int)
			INPLACE_REAL_OPERATOR(/=,double)
			INPLACE_REAL_OPERATOR(/=,value_type)
			INPLACE_COMPLEX_MULT(int)
			INPLACE_COMPLEX_MULT(double)
			INPLACE_COMPLEX_DIV(int)
			INPLACE_COMPLEX_DIV(double)
		private:
			template <class Complex>
			complex &divide_by_complex(const Complex &other)
			{
				if (other.real() == 0 && other.imag() == 0) {
					piranha_throw(zero_division_error,"cannot divide by zero");
				}
				// This is the divisor, i.e. the square of absolute value of other.
				value_type div(other.real());
				div *= other.real();
				div += other.imag() * other.imag();
				// The numerator looks like a multiplication with opposite signs.
				const value_type tmp1 = m_imag * other.imag(), tmp2 = m_real * other.imag();
				m_imag *= other.real();
				m_imag -= tmp2;
				m_real *= other.real();
				m_real += tmp1;
				// Now divide by divisor.
				m_real /= div;
				m_imag /= div;
				return *this;
			}
			template <class Complex>
			complex &mult_by_complex(const Complex &other)
			{
				const value_type tmp1(m_imag * other.imag()), tmp2(m_real * other.imag());
				// NOTE: we do imag first because if we modify real now, then it screws up the computation.
				// m_imag is not used anymore as rhs from this point onwards.
				m_imag *= other.real();
				m_imag += tmp2;
				m_real *= other.real();
				m_real -= tmp1;
				return *this;
			}
			void construct_from_string(const char *str)
			{
				string tmp(str);
				// First let's trim the input string.
				boost::trim(tmp);
				// For string to be a valid complex number, it must consist of a least 5 chars
				// (2 brackets, a comma and the two numbers).
				if (tmp.size() < 5) {
					piranha_throw(value_error,"invalid string input");
				}
				// Next we split the input string into two parts, separated by the comma.
				vector<string> split_v;
				boost::split(split_v,tmp,boost::is_any_of(","));
				// Now we check that the input string is well formed, i.e., it contains one comma,
				// it starts and ends with round brackets and each element of the split vector contains at
				// least two characters.
				if (split_v.size() != 2 || split_v[0].size() < 2 || split_v[1].size() < 2 ||
					split_v[0][0] != '(' || split_v[1][split_v[1].size() - 1] != ')') {
					piranha_throw(value_error,"invalid string input");
				}
				// Finally, we try to build the two components from the split vector,
				// discarding the positions in which there are the brackets.
				m_real = value_type(string(&split_v[0][1],&split_v[0][split_v[0].size()]));
				m_imag = value_type(string(&split_v[1][0],&split_v[1][split_v[1].size() - 1]));
			}
			value_type	m_real;
			value_type	m_imag;
	};

	#undef INPLACE_REAL_OPERATOR
	#undef INPLACE_COMPLEX_OPERATOR
	#undef INPLACE_COMPLEX_MULT
	#undef INPLACE_COMPLEX_DIV
	#undef REAL_EQUALITY
	#undef COMPLEX_EQUALITY

	/// Overload of standard swap function for std::complex<piranha::mp_rational>.
	/**
	 * Will use the swap() method internally.
	 */
	inline void swap(complex<piranha::mp_rational> &qc1, complex<piranha::mp_rational> &qc2)
	{
		qc1.swap(qc2);
	}

	/// Overload out stream operator<< for std::complex<piranha::mp_rational>.
	inline ostream &operator<<(ostream &o, const complex<piranha::mp_rational> &qc)
	{
		o << '(' << qc.m_real << ',' << qc.m_imag << ')';
		return o;
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
}

#endif
