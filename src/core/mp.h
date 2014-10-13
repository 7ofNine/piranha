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

#ifndef PIRANHA_MP_H
#define PIRANHA_MP_H

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>

#include "math.h"
// For now we support only GMP.
#include "mp/piranha_gmp.h"

// Function overloads that go into the piranha namespace.
namespace piranha
{
	/* RATIONAL CLASS OVERLOADS */
	/* ------------------------ */

	/// Overload of out stream operator<< for piranha::mp_rational.
	inline std::ostream &operator<<(std::ostream &o, const mp_rational &q)
	{
		return q.print(o);
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
	inline std::size_t hash_value(const mp_rational &q)
	{
                return q.hash();
	}

	/// Overload root function for piranha::mp_rational.
	inline mp_rational root(const int &n, const mp_rational &q)
	{
		return q.root(n);
	}

	/// Overload r_factorial for piranha::mp_rational.
	inline mp_rational r_factorial(const mp_rational &q, const int &n)
	{
		return q.r_factorial(n);
	}

	/// Overload f_factorial for piranha::mp_rational.
	inline mp_rational f_factorial(const mp_rational &q, const int &n)
	{
		return q.f_factorial(n);
	}

	/// Overload binomial coefficient (choose) function for piranha::mp_rational.
	/**
	 * @see piranha::mp_rational::choose.
	 */
	inline mp_rational choose(const mp_rational &q, const int &k)
	{
		return q.choose(k);
	}

	/// Overload in stream operator>> for std::complex<piranha::mp_rational>.
	inline std::istream &operator>>(std::istream &i, std::complex<mp_rational> &qc)
	{
		std::string tmp_str;
		std::getline(i,tmp_str);
		std::complex<mp_rational> tmp(tmp_str);
		qc.swap(tmp);
		return i;
	}

	/// Overload hash_value function for std::complex<piranha::mp_rational>.
	inline std::size_t hash_value(const std::complex<mp_rational> &qc)
	{
		return qc.hash();
	}

	/// Overload root function for std::complex<piranha::mp_rational>.
	inline std::complex<mp_rational> root(const int &n, const std::complex<mp_rational> &qc)
	{
		return qc.root(n);
	}

	/// Overload r_factorial for std::complex<piranha::mp_rational>.
	inline std::complex<mp_rational> r_factorial(const std::complex<mp_rational> &c, const int &n)
	{
		return c.r_factorial(n);
	}

	/// Overload f_factorial for std::complex<piranha::mp_rational>.
	inline std::complex<mp_rational> f_factorial(const std::complex<mp_rational> &c, const int &n)
	{
		return c.f_factorial(n);
	}

	/// Overload binomial coefficient (choose) function for std::complex<piranha::mp_rational>.
	/**
	 * @see std::complex<piranha::mp_rational>::choose.
	 */
	inline std::complex<mp_rational> choose(const std::complex<mp_rational> &c, const int &k)
	{
		return c.choose(k);
	}

	/// Check if input argument can be represented exactly as an int type.
	inline bool is_integer(const mp_rational &q)
	{
		try {
			int n = q.to_int();
			(void)n;
			return true;
		} catch (const value_error &) {
			// This means that x is not an integer.
		} catch (const std::overflow_error &) {
			// This means that x is an integer but overflows int range.
		}
		return false;
	}

	/// Bessel function of the first kind.
	inline mp_rational besselJ(const int &order, const mp_rational &arg)
	{
		if (arg != 0) {
				PIRANHA_THROW(value_error,"cannot compute Bessel function of non-zero rational");
		}
		return (order == 0) ? mp_rational(1) : mp_rational(0);
	}


	/* INTEGER CLASS OVERLOADS */
	/* ----------------------- */

	/// Overload of out stream operator<< for piranha::mp_integer.
	inline std::ostream &operator<<(std::ostream &o, const mp_integer &z)
	{
		return z.print(o);
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
	inline std::size_t hash_value(const mp_integer &z)
	{
                return z.hash();
	}

	/// Overload root function for piranha::mp_integer.
	inline mp_integer root(const int &n, const mp_integer &z)
	{
		return z.root(n);
	}

	/// Overload factorial function for piranha::mp_integer.
	/**
	 * @see piranha::mp_integer::factorial.
	 */
	inline mp_integer factorial(const mp_integer &z)
	{
		return z.factorial();
	}

	/// Overload double factorial function for piranha::mp_integer.
	/**
	 * @see piranha::mp_integer::double_factorial.
	 */
	inline mp_integer double_factorial(const mp_integer &z)
	{
		return z.double_factorial();
	}

	/// Overload r_factorial for piranha::mp_integer.
	inline mp_integer r_factorial(const mp_integer &z, const int &n)
	{
		return z.r_factorial(n);
	}

	/// Overload f_factorial for piranha::mp_integer.
	inline mp_integer f_factorial(const mp_integer &z, const int &n)
	{
		return z.f_factorial(n);
	}

	/// Overload binomial coefficient (choose) function for piranha::mp_integer.
	/**
	 * @see piranha::mp_integer::choose.
	 */
	inline mp_integer choose(const mp_integer &n, const int &k)
	{
		return n.choose(k);
	}

	/// Multiply-accumulate overload for piranha::mp_integer.
	/**
	 * Will call mp_integer::multiply_accumulate().
	 */
	inline void multiply_accumulate(mp_integer &x, const mp_integer &y, const mp_integer &z)
	{
		x.multiply_accumulate(y,z);
	}

	/// Overload in stream operator>> for std::complex<piranha::mp_integer>.
	inline std::istream &operator>>(std::istream &i, std::complex<mp_integer> &zc)
	{
		std::string tmp_str;
		std::getline(i,tmp_str);
		std::complex<mp_integer> tmp(tmp_str);
		zc.swap(tmp);
		return i;
	}

	/// Overload hash_value function for std::complex<piranha::mp_integer>.
	inline std::size_t hash_value(const std::complex<mp_integer> &zc)
	{
		return zc.hash();
	}

	/// Overload root function for std::complex<piranha::mp_integer>.
	inline std::complex<mp_integer> root(const int &n, const std::complex<mp_integer> &zc)
	{
		return zc.root(n);
	}

	/// Overload r_factorial for std::complex<piranha::mp_integer>.
	inline std::complex<mp_integer> r_factorial(const std::complex<mp_integer> &c, const int &n)
	{
		return c.r_factorial(n);
	}

	/// Overload f_factorial for std::complex<piranha::mp_integer>.
	inline std::complex<mp_integer> f_factorial(const std::complex<mp_integer> &c, const int &n)
	{
		return c.f_factorial(n);
	}

	/// Overload binomial coefficient (choose) function for std::complex<piranha::mp_integer>.
	/**
	 * @see std::complex<piranha::mp_integer>::choose.
	 */
	inline std::complex<mp_integer> choose(const std::complex<mp_integer> &c, const int &k)
	{
		return c.choose(k);
	}

	/// Check if input argument can be represented exactly as an int type.
	inline bool is_integer(const mp_integer &z)
	{
		try {
			int n = z.to_int();
			(void)n;
			return true;
		} catch (const std::overflow_error &) {
			// This means that x is an integer but overflows int range.
		}
		return false;
	}

	/// Condon-Shortley phase.
	inline int cs_phase(const mp_integer &z)
	{
		// NOTE: more efficient implementation here, probably with GMP functions?
		if (z % 2 != 0) {
			return -1;
		} else {
			return 1;
		}
	}

	/// Bessel function of the first kind.
	inline mp_integer besselJ(const int &order, const mp_integer &arg)
	{
		if (arg != 0) {
				PIRANHA_THROW(value_error,"cannot compute Bessel function of non-zero integer");
		}
		return (order == 0) ? mp_integer(1) : mp_integer(0);
	}
}

namespace std
{
	/* RATIONAL CLASS OVERLOADS */
	/* ------------------------ */

	/// Overload standard swap function for piranha::mp_rational.
	/**
	 * Will use piranha::mp_rational::swap() internally.
	 * @see piranha::mp_rational::swap().
	 */
	inline void swap(piranha::mp_rational &q1, piranha::mp_rational &q2)
	{
		q1.swap(q2);
	}

	#define STD_POW_OVERLOAD(ret_type, expo_type) \
	/** \brief Overload standard power function for ret_type argument and expo_type exponent. */ \
	inline ret_type pow(const ret_type &x, const expo_type &e) \
	{ \
		return x.pow(e); \
	}

	STD_POW_OVERLOAD(piranha::mp_rational,double)
	STD_POW_OVERLOAD(piranha::mp_rational,piranha::mp_rational)

	/// Overload standard abs function for piranha::mp_rational.
	/**
	 * @see piranha::mp_rational::abs.
	 */
	inline piranha::mp_rational abs(const piranha::mp_rational &q)
	{
		return q.abs();
	}

	/// Overload of standard swap function for std::complex<piranha::mp_rational>.
	/**
	 * Will use the swap() method internally.
	 */
	inline void swap(complex<piranha::mp_rational> &qc1, complex<piranha::mp_rational> &qc2)
	{
		qc1.swap(qc2);
	}

	STD_POW_OVERLOAD(std::complex<piranha::mp_rational>,double)
	STD_POW_OVERLOAD(std::complex<piranha::mp_rational>,piranha::mp_rational)

	/// Overload standard abs function for std::complex<piranha::mp_rational>.
	/**
	 * @see complex_generic_mp_container::abs()
	 */
	inline piranha::mp_rational abs(const complex<piranha::mp_rational> &qc)
	{
		return qc.abs();
	}

	/* INTEGER CLASS OVERLOADS */
	/* ------------------------ */

	/// Overload standard swap function for piranha::mp_integer.
	/**
	 * Will use piranha::mp_integer::swap() internally.
	 * @see piranha::mp_integer::swap().
	 */
	inline void swap(piranha::mp_integer &z1, piranha::mp_integer &z2)
	{
		z1.swap(z2);
	}

	STD_POW_OVERLOAD(piranha::mp_integer,double)
	STD_POW_OVERLOAD(piranha::mp_integer,piranha::mp_rational)

	/// Overload standard abs function for piranha::mp_integer.
	/**
	 * @see piranha::mp_integer::abs.
	 */
	inline piranha::mp_integer abs(const piranha::mp_integer &z)
	{
		return z.abs();
	}

	/// Overload of standard swap function for std::complex<piranha::mp_integer>.
	/**
	 * Will use the swap() method internally.
	 */
	inline void swap(complex<piranha::mp_integer> &zc1, complex<piranha::mp_integer> &zc2)
	{
		zc1.swap(zc2);
	}

	STD_POW_OVERLOAD(std::complex<piranha::mp_integer>,double)
	STD_POW_OVERLOAD(std::complex<piranha::mp_integer>,piranha::mp_rational)

	/// Overload standard abs function for std::complex<piranha::mp_integer>.
	/**
	 * @see complex_generic_mp_container::abs()
	 */
	inline piranha::mp_integer abs(const complex<piranha::mp_integer> &zc)
	{
		return zc.abs();
	}

	/* OVERLOADS OF POD NUMERICAL TYPES MATHS */
	/* -------------------------------------- */

	/// Overload std::pow for double and piranha::mp_rational arguments.
	inline double pow(const double &x, const piranha::mp_rational &q)
	{
		return pow(x,q.to_double());
	}

	/// Overload std::pow for std::complex<double> and piranha::mp_rational arguments.
	inline complex<double> pow(const complex<double> &c, const piranha::mp_rational &q)
	{
		return pow(c,q.to_double());
	}

	#undef STD_POW_OVERLOAD
}

#endif
