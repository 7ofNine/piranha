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
#include <boost/functional/hash.hpp>
#include <boost/operators.hpp>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <string>

#include "../exceptions.h"

namespace piranha
{
	/// Multiprecision rational class.
	/**
	 * Wraps a GMP mpq_class. Interoperability with C++ int and double types is provided through constructors and
	 * in-place mathematical operators against arbitrary types that provide forwarding to the underlying
	 * mpq_class type. Additional operators are provided through inheritance from the Boost operator library.	
	 */
	// TODO: catch std::invalid_argument when initialising from invalid string
	// and rethrow as value_error.
	class mp_rational:
		boost::ordered_field_operators<mp_rational,
		boost::ordered_field_operators<mp_rational, int,
		boost::ordered_field_operators<mp_rational, double
		> > >
	{
			friend std::ostream &operator<<(std::ostream &, const mp_rational &);
			friend std::istream &operator>>(std::istream &, mp_rational &);
		public:
			/// Default constructor.
			/**
			 * Initialises value to zero.
			 */
			explicit mp_rational():m_value(0) {}
			/// Constructor from std::string.
			/**
			 * Will raise a piranha::value_error exception if string is not valid. Valid strings include
			 * integer numbers and rationals, even with negative signs (e.g., "3", "-10", "3/4", "5/-6", etc.).
			 */
			explicit mp_rational(const std::string &str):m_value(str.c_str())
			{
				canonicalize();
			}
			/// Constructor from C string.
			/**
			 * @see mp_rational(const std::string &).
			 */
			explicit mp_rational(const char *str):m_value(str)
			{
				canonicalize();
			}
			/// Constructor from integer.
			explicit mp_rational(const int &n):m_value(n) {}
			/// Constructor from integer numerator and denominator.
			explicit mp_rational(const int &n, const int &d):m_value(n,d)
			{
				canonicalize();
			}
			/// Constructor from double.
			explicit mp_rational(const double &x):m_value(x) {}
			/// Copy constructor.
			mp_rational(const mp_rational &other):m_value(other.m_value) {}
			/// Assignment operator.
			mp_rational &operator=(const mp_rational &other)
			{
				m_value = other.m_value;
				return *this;
			}
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
			mp_rational &operator/=(const mp_rational &other)
			{
				if (other.m_value == 0) {
					throw division_by_zero();
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
			/// Assignment operator for arbitrary type.
			template <class T>
			mp_rational &operator=(const T &other)
			{
				m_value = other;
				return *this;
			}
			/// Equality operator with arbitrary type.
			template <class T>
			bool operator==(const T &other) const
			{
				return (m_value == other);
			}
			/// Comparison operator with arbitrary type.
			template <class T>
			bool operator<(const T &other) const
			{
				return (m_value < other);
			}
			/// In-place addition with arbitrary type.
			template <class T>
			mp_rational &operator+=(const T &other)
			{
				m_value += other;
				return *this;
			}
			/// In-place subtraction with arbitrary type.
			template <class T>
			mp_rational &operator-=(const T &other)
			{
				m_value -= other;
				return *this;
			}
			/// In-place multiplication with arbitrary type.
			template <class T>
			mp_rational &operator*=(const T &other)
			{
				m_value *= other;
				return *this;
			}
			/// In-place division with arbitrary type.
			template <class T>
			mp_rational &operator/=(const T &other)
			{
				if (other == 0) {
					throw division_by_zero();
				}
				m_value /= other;
				return *this;
			}
		private:
			void canonicalize()
			{
				mpq_canonicalize(m_value.get_mpq_t());
			}
			mpq_class m_value;
	};

	/// Overload out stream operator<< for piranha::mp_rational.
	inline std::ostream &operator<<(std::ostream &o, const mp_rational &q)
	{
		o << q.m_value;
		return o;
	}

	/// Overload in stream operator>> for piranha::mp_rational.
	inline std::istream &operator>>(std::istream &i, mp_rational &q)
	{
		i >> q.m_value;
		return i;
	}

	/// Overload hash_value function for use with boost::hash.
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
}

#endif
