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

#ifndef PIRANHA_GENERIC_MP_CONTAINER_H
#define PIRANHA_GENERIC_MP_CONTAINER_H

#include <boost/algorithm/string.hpp>
#include <boost/operators.hpp>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	// Macros for operators.
	#define ASSIGN_METHOD(type) \
	/** \brief Assign type. Will call operator= internally. */ \
	Derived &assign(const type &x) \
	{ \
		return this->operator=(x); \
	}
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
	#define INPLACE_OPERATOR(op,type) \
	/** \brief In-place operator op for type. */ \
	Derived &operator op(const type &other) \
	{ \
		m_value op other; \
		return *derived_cast; \
	}
	#define INPLACE_DIVISION(type) \
	/** \brief In-place division with type. */ \
	Derived &operator/=(const type &other) \
	{ \
		if (other == 0) { \
			piranha_throw(zero_division_error,"cannot divide by zero"); \
		} \
		m_value /= other; \
		return *derived_cast; \
	}

	/// Generic container helper for multi-precision classes.
	/**
	 * Can be used to wrap around GMP C++ classes. It provides basic mathematical operators
	 * against int, double and the Derived class. Uses the CRTP internally.
	 */
	template <class T, class Derived>
	class generic_mp_container:
		boost::ordered_field_operators<Derived,
		boost::ordered_field_operators<Derived, int,
		boost::ordered_field_operators<Derived, double
		> > >
	{
			template <class T2, class Derived2>
			friend std::ostream &operator<<(std::ostream &, const generic_mp_container<T2,Derived2> &);
		public:
			/// Default constructor.
			/**
			 * Will initialise content to zero.
			 */
			explicit generic_mp_container():m_value(0) {}
			/// Constructor from integer.
			explicit generic_mp_container(const int &n):m_value(n) {}
			/// Constructor from double.
			explicit generic_mp_container(const double &x):m_value(x) {}
			/// Constructor from string.
			explicit generic_mp_container(const std::string &str):m_value(str) {}
			/// Constructor from C string.
			explicit generic_mp_container(const char *str):m_value(str) {}
			/// Equality operator.
			bool operator==(const Derived &other) const
			{
				return (m_value == other.m_value);
			}
			/// Comparison operator.
			bool operator<(const Derived &other) const
			{
				return (m_value < other.m_value);
			}
			/// In-place addition.
			Derived &operator+=(const Derived &other)
			{
				m_value += other.m_value;
				return *derived_cast;
			}
			/// In-place subtraction.
			Derived &operator-=(const Derived &other)
			{
				m_value -= other.m_value;
				return *derived_cast;
			}
			/// In-place multiplication.
			Derived &operator*=(const Derived &other)
			{
				m_value *= other.m_value;
				return *derived_cast;
			}
			/// In-place division.
			/**
			 * @throws zero_division_error if dividing by zero.
			 */
			Derived &operator/=(const Derived &other)
			{
				if (other.m_value == 0) {
					piranha_throw(zero_division_error,"cannot divide by zero");
				}
				m_value /= other.m_value;
				return *derived_cast;
			}
			// Operators against plain old numerical types.
			EQUALITY_OPERATOR(int)
			EQUALITY_OPERATOR(double)
			COMPARISON_OPERATOR(int)
			COMPARISON_OPERATOR(double)
			INPLACE_OPERATOR(=,double)
			INPLACE_OPERATOR(=,int)
			INPLACE_OPERATOR(+=,int)
			INPLACE_OPERATOR(+=,double)
			INPLACE_OPERATOR(-=,int)
			INPLACE_OPERATOR(-=,double)
			INPLACE_OPERATOR(*=,int)
			INPLACE_OPERATOR(*=,double)
			INPLACE_DIVISION(int)
			INPLACE_DIVISION(double)
		protected:
			// These are provided because the derived class will not inherit the assignment operator.
			// By using this method it is possible to call the assignment operators here defined
			// in the derived class.
			ASSIGN_METHOD(int)
			ASSIGN_METHOD(double)
		protected:
			T m_value;
	};

	#undef EQUALITY_OPERATOR
	#undef COMPARISON_OPERATOR
	#undef INPLACE_OPERATOR
	#undef INPLACE_DIVISION

	/// Overload of out stream operator<<.
	template <class T, class Derived>
	inline std::ostream &operator<<(std::ostream &o, const generic_mp_container<T,Derived> &c)
	{
		o << c.m_value;
		return o;
	}

	// Useful macros for operators.
	#define INPLACE_REAL_OPERATOR(op,type) \
	/** \brief In-place operator op for type. */ \
	Derived &operator op(const type &x) \
	{ \
		m_real op x; \
		return *derived_cast; \
	}
	#define INPLACE_REAL_OPERATOR2(op,type) \
	/** \brief In-place operator op for type. */ \
	Derived &operator op(const type &x) \
	{ \
		m_real op x; \
		m_imag op x; \
		return *derived_cast; \
	}
	#define INPLACE_REAL_ASSIGN(type) \
	/** \brief Assign type. */ \
	Derived &operator=(const type &x) \
	{ \
		m_real = x; \
		m_imag = 0; \
		return *derived_cast; \
	}
	#define INPLACE_COMPLEX_OPERATOR(op,type) \
	/** \brief In-place operator op for complex type. */ \
	Derived &operator op(const std::complex< type > &c) \
	{ \
		m_real op c.real(); \
		m_imag op c.imag(); \
		return *derived_cast; \
	}
	#define INPLACE_COMPLEX_MULT(type) \
	/** \brief In-place operator *= for complex type. */ \
	Derived &operator*=(const std::complex< type > &c) \
	{ \
		return mult_by_complex(c); \
	}
	#define INPLACE_COMPLEX_DIV(type) \
	/** \brief In-place operator /= for complex type. */ \
	Derived &operator/=(const std::complex< type > &c) \
	{ \
		return divide_by_complex(c); \
	}
	#define REAL_EQUALITY(type) \
	/** \brief Equality operator for type. */ \
	bool operator==(const type &x) const \
	{ \
		return (m_real == x && m_imag == 0); \
	}
	#define COMPLEX_EQUALITY(type) \
	/** \brief Equality operator for complex type. */ \
	bool operator==(const std::complex< type > &c) const \
	{ \
		return (m_real == c.real() && m_imag == c.imag()); \
	}

	/// Generic container helper for complex multi-precision classes.
	/**
	 * Can be used to help building the complex counterparts of wrapped GMP classes.
	 * Provides basic mathematical interoperability with int, double, complex int, complex
	 * double and the real counterpart. The constructors will forward the construction to the
	 * underlying real class. Uses the CRTP.
	 */
	template <class T, class Derived>
	class complex_generic_mp_container:
		boost::field_operators<Derived,
		boost::field_operators<Derived, int,
		boost::field_operators<Derived, double,
		boost::field_operators<Derived, T,
		boost::field_operators<Derived, std::complex<int>,
		boost::field_operators<Derived, std::complex<double>
		> > > > > >
	{
		public:
			/// Default constructor.
			explicit complex_generic_mp_container(): m_real(0),m_imag(0) {}
			/// Constructor from integer.
			explicit complex_generic_mp_container(const int &n): m_real(n),m_imag(0) {}
			/// Constructor from integer real and imaginary parts..
			explicit complex_generic_mp_container(const int &r, const int &i): m_real(r),m_imag(i) {}
			/// Constructor from complex integer.
			explicit complex_generic_mp_container(const std::complex<int> &c):
				m_real(c.real()),m_imag(c.imag()) {}
			/// Constructor from double.
			explicit complex_generic_mp_container(const double &x): m_real(x),m_imag(0) {}
			/// Constructor from double real and imaginary parts.
			explicit complex_generic_mp_container(const double &r, const double &i): m_real(r),m_imag(i) {}
			/// Constructor from complex double.
			explicit complex_generic_mp_container(const std::complex<double> &c):
				m_real(c.real()),m_imag(c.imag()) {}
			/// Constructor from T.
			explicit complex_generic_mp_container(const T &x): m_real(x),m_imag(0) {}
			/// Constructor from real and imaginary parts of type value_type.
			explicit complex_generic_mp_container(const T &r, const T &i): m_real(r),m_imag(i) {}
			/// Constructor from std::string.
			/**
			 * Will raise a value_error exception if string is not valid. A valid string is of the form
			 * "(real,imag)", where "real" and "imag" are strings that can construct successfully an
			 * instance of type T.
			 * @throws value_error if string is invalid.
			 */
			explicit complex_generic_mp_container(const std::string &str): m_real(),m_imag()
			{
				construct_from_string(str.c_str());
			}
			/// Constructor from C string.
			/**
			 * @see complex_generic_mp_container(const std::string &).
			 */
			explicit complex_generic_mp_container(const char *str): m_real(),m_imag()
			{
				construct_from_string(str);
			}
			/// Get const reference to the real part.
			const T &real() const
			{
				return m_real;
			}
			/// Get const reference to the imaginary part.
			const T &imag() const
			{
				return m_imag;
			}
			/// Equality operator.
			bool operator==(const Derived &other) const
			{
				return (m_real == other.m_real && m_imag == other.m_imag);
			}
			// Math operators.
			/// In-place addition.
			Derived &operator+=(const Derived &other)
			{
				m_real += other.m_real;
				m_imag += other.m_imag;
				return *derived_cast;
			}
			/// In-place subtraction.
			Derived &operator-=(const Derived &other)
			{
				m_real -= other.m_real;
				m_imag -= other.m_imag;
				return *derived_cast;
			}
			/// In-place multiplication.
			Derived &operator*=(const Derived &other)
			{
				return mult_by_complex(other);
			}
			/// In-place division.
			Derived &operator/=(const Derived &other)
			{
				return divide_by_complex(other);
			}
			// Maths for other types.
			REAL_EQUALITY(int)
			REAL_EQUALITY(double)
			REAL_EQUALITY(T)
			COMPLEX_EQUALITY(int)
			COMPLEX_EQUALITY(double)
			INPLACE_REAL_ASSIGN(int)
			INPLACE_REAL_ASSIGN(double)
			INPLACE_REAL_ASSIGN(T)
			INPLACE_COMPLEX_OPERATOR(=,int)
			INPLACE_COMPLEX_OPERATOR(=,double)
			INPLACE_REAL_OPERATOR(+=,int)
			INPLACE_REAL_OPERATOR(+=,double)
			INPLACE_REAL_OPERATOR(+=,T)
			INPLACE_COMPLEX_OPERATOR(+=,int)
			INPLACE_COMPLEX_OPERATOR(+=,double)
			INPLACE_REAL_OPERATOR(-=,int)
			INPLACE_REAL_OPERATOR(-=,double)
			INPLACE_REAL_OPERATOR(-=,T)
			INPLACE_COMPLEX_OPERATOR(-=,int)
			INPLACE_COMPLEX_OPERATOR(-=,double)
			INPLACE_REAL_OPERATOR2(*=,int)
			INPLACE_REAL_OPERATOR2(*=,double)
			INPLACE_REAL_OPERATOR2(*=,T)
			INPLACE_REAL_OPERATOR2(/=,int)
			INPLACE_REAL_OPERATOR2(/=,double)
			INPLACE_REAL_OPERATOR2(/=,T)
			INPLACE_COMPLEX_MULT(int)
			INPLACE_COMPLEX_MULT(double)
			INPLACE_COMPLEX_DIV(int)
			INPLACE_COMPLEX_DIV(double)
			/// Mathematical inversion.
			Derived invert() const
			{
				Derived retval(*derived_const_cast);
				retval.m_imag.negate();
				const T div = abs2();
				retval.m_real /= div;
				retval.m_imag /= div;
				return retval;
			}
		private:
			// Square of absolute value.
			T abs2() const
			{
				// NOTE: rewrite this in terms of multadd, when implemented.
				T retval(m_real);
				retval *= m_real;
				retval += m_imag * m_imag;
				return retval;
			}
			template <class Complex>
			Derived &divide_by_complex(const Complex &other)
			{
				// NOTE: rewrite this in terms of multadd, when implemented.
				if (other.real() == 0 && other.imag() == 0) {
					piranha_throw(zero_division_error,"cannot divide by zero");
				}
				// This is the divisor, i.e. the square of absolute value of other.
				T div(other.real());
				div *= other.real();
				div += other.imag() * other.imag();
				// The numerator looks like a multiplication with opposite signs.
				const T tmp1 = m_imag * other.imag(), tmp2 = m_real * other.imag();
				m_imag *= other.real();
				m_imag -= tmp2;
				m_real *= other.real();
				m_real += tmp1;
				// Now divide by divisor.
				m_real /= div;
				m_imag /= div;
				return *derived_cast;
			}
			template <class Complex>
			Derived &mult_by_complex(const Complex &other)
			{
				// NOTE: rewrite this in terms of multadd, when implemented?
				const T tmp1(m_imag * other.imag()), tmp2(m_real * other.imag());
				// NOTE: we do imag first because if we modify real now, then it screws up the computation.
				// m_imag is not used anymore as rhs from this point onwards.
				m_imag *= other.real();
				m_imag += tmp2;
				m_real *= other.real();
				m_real -= tmp1;
				return *derived_cast;
			}
			void construct_from_string(const char *str)
			{
				std::string tmp(str);
				// First let's trim the input string.
				boost::trim(tmp);
				// Next we split the input string into two parts, separated by the comma.
				std::vector<std::string> split_v;
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
				m_real = T(std::string(&split_v[0][1],&split_v[0][split_v[0].size()]));
				m_imag = T(std::string(&split_v[1][0],&split_v[1][split_v[1].size() - 1]));
			}
		protected:
			ASSIGN_METHOD(T)
			ASSIGN_METHOD(int)
			ASSIGN_METHOD(double)
			ASSIGN_METHOD(std::complex<int>)
			ASSIGN_METHOD(std::complex<double>)
		protected:
			T	m_real;
			T	m_imag;
	};

	#undef INPLACE_REAL_OPERATOR
	#undef INPLACE_COMPLEX_OPERATOR
	#undef INPLACE_COMPLEX_MULT
	#undef INPLACE_COMPLEX_DIV
	#undef REAL_EQUALITY
	#undef COMPLEX_EQUALITY
	#undef ASSIGN_METHOD

	/// Overload out stream operator<<
	template <class T, class Derived>
	inline std::ostream &operator<<(std::ostream &o, const complex_generic_mp_container<T,Derived> &c)
	{
		o << '(' << c.real() << ',' << c.imag() << ')';
		return o;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
