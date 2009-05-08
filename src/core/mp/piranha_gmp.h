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
#include "generic_mp_container.h"

namespace piranha
{
	#define ASSIGNMENT_OPERATOR(self, type) \
	/** \brief Assignment operator for type. */ \
	self &operator=(const type &x) \
	{ \
		return this->assign(x); \
	}

	/// Multiprecision rational class.
	/**
	 * Wraps a GMP mpq_class.
	 */
	class mp_rational: public generic_mp_container<mpq_class,mp_rational>
	{
			typedef generic_mp_container<mpq_class,mp_rational> ancestor;
		public:
			/// Default constructor.
			/**
			 * Initialises value to zero.
			 */
			explicit mp_rational():ancestor() {}
			/// Constructor from std::string.
			/**
			 * Will raise a value_error exception if string is not valid. Valid strings include
			 * integer numbers and rationals, even with negative signs (e.g., "3", "-10", "3/4", "5/-6", etc.).
			 * @throws value_error if string is invalid.
			 * @throws zero_division_error if string is valid but denominator is zero.
			 */
			explicit mp_rational(const std::string &str):ancestor()
			{
				construct_from_string(str.c_str());
			}
			/// Constructor from C string.
			/**
			 * @see mp_rational(const std::string &).
			 */
			explicit mp_rational(const char *str):ancestor()
			{
				construct_from_string(str);
			}
			/// Constructor from integer.
			explicit mp_rational(const int &n):ancestor(n) {}
			/// Constructor from integer numerator and denominator.
			/**
			 * @throws zero_division_error if denominator is zero.
			 */
			explicit mp_rational(const int &n, const int &d):ancestor()
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
			explicit mp_rational(const double &x):ancestor(x) {}
			ASSIGNMENT_OPERATOR(mp_rational,int)
			ASSIGNMENT_OPERATOR(mp_rational,double)
			/// Cast to double.
			/**
			 * Uses the get_d() GMP routine internally.
			 */
			operator double() const
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
			operator int() const
			{
				if (m_value.get_den() != 1) {
					piranha_throw(value_error,"cannot convert rational to integer if denominator is non-unitary");
				}
				if (!m_value.get_num().fits_sint_p()) {
					piranha_throw(std::overflow_error,"numerator is too large while converting rational to integer");
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
	};

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

	/// Complex counterpart of piranha::mp_rational.
	/**
	 * Stores two piranha::mp_rational members as real and imaginary part.
	 */
	template <>
	class complex<piranha::mp_rational>:
		public piranha::complex_generic_mp_container<piranha::mp_rational,complex<piranha::mp_rational> >
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
			ASSIGNMENT_OPERATOR(complex,int)
			ASSIGNMENT_OPERATOR(complex,double)
			ASSIGNMENT_OPERATOR(complex,complex<int>)
			ASSIGNMENT_OPERATOR(complex,complex<double>)
			ASSIGNMENT_OPERATOR(complex,value_type)
			/// Cast to complex double.
			/**
			 * @see piranha::mp_rational::operator double.
			 */
			operator complex<double>() const
			{
				return complex<double>((double)m_real,(double)m_imag);
			}
			/// Cast to complex int.
			/**
			 * @see piranha::mp_rational::operator int.
			 */
			operator complex<int>() const
			{
				return complex<int>((int)m_real,(int)m_imag);
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
			complex pow_int(const int &n) const {
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
			complex pow_double(const double &y) const {
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

	#undef ASSIGNMENT_OPERATOR
}

#endif
