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

#ifndef PIRANHA_MATH_H
#define PIRANHA_MATH_H

#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/static_assert.hpp>
#include <cmath>
#include <complex>

#include "config.h"
#include "exceptions.h"

namespace piranha
{
	/// Meta-programmed functor for the calculation of base-2 logarithm.
	/**
	 * Result is retrieved through the lg::value const static member.
	 */
	template <int N>
	struct lg {
		p_static_check(N > 0 && (N % 2) == 0, "N is not a power of 2.");
		/// Value of the base-2 logarithm of N.
		static const size_t value = lg < (N >> 1) >::value + 1;
	};

	template <>
	struct lg<1> {
		static const size_t value = 0;
	};

	// Integer base-2 logarithm. Like lg, but won't error out if N is not a power of 2.
	template <int N>
	struct ilg {
		p_static_check(N > 0, "N must be positive.");
		static const size_t value = ilg < (N >> 1) >::value + 1;
	};

	template <>
	struct ilg<1> {
		static const size_t value = 0;
	};

	template <long int Xk, long int N, bool Stop, long int Xkm1 = 0>
	struct isqrt_impl {
		// We should be safe here, since N == Xk initially and then Xk will
		// decrease at every iteration for every meaningful value of cache size.
		static const long int Xkp1 = (Xk + N/Xk) / 2;
		static const long int value = isqrt_impl<Xkp1,N, Xk - Xkp1 <= 0, Xk>::value;
	};

	template <long int Xk, long int N, long int Xkm1>
	struct isqrt_impl<Xk, N, true, Xkm1> {
		static const long int value = (Xkm1 < Xk) ? Xkm1 : Xk;
	};

	// Integer square root. See http://en.wikipedia.org/wiki/Integer_square_root for implementation
	// (Newton's method).
	template <long int N>
	struct isqrt {
		p_static_check(N > 0, "Invalid value for isqrt.");
		static const long int value_ = isqrt_impl<N,N,false>::value;
		p_static_check(value_ >= 0, "Invalid result for isqrt.");
		static const size_t value = static_cast<size_t>(value_);
	};

	template <>
	struct isqrt<0> {
		static const size_t value = 0;
	};

	template <class T>
	inline int sign(const T &x)
	{
		if (x >= 0) {
			return 1;
		} else {
			return -1;
		}
	}

	/// Condon-Shortley phase.
	inline int cs_phase(const int &n)
	{
		if (n & 1) {
			return -1;
		} else {
			return 1;
		}
	}

	template <class T>
	static inline void factorial_check(const T &x)
	{
		if (x < 0) {
			piranha_throw(value_error,
				"please use a non-negative integer as argument for factorials");
		}
	}

	/// Factorial.
	/**
	 * @throws value_error if i is negative.
	 * @throws std::overflow_error if result is too big to be represented in double-precision.
	 */
	inline double factorial(const int &i)
	{
		factorial_check(i);
		return boost::math::factorial<double>(static_cast<unsigned>(i));
	}

	/// Double factorial.
	/**
	 * @throws value_error if i is negative.
	 * @throws std::overflow_error if result is too big to be represented in double-precision.
	 */
	inline double double_factorial(const int &i)
	{
		factorial_check(i);
		return boost::math::double_factorial<double>(static_cast<unsigned>(i));
	}

	template <class Number>
	inline Number generic_r_factorial(const Number &x, const int &n)
	{
		factorial_check(n);
		Number retval(1);
		for (int i = 0; i < n; ++i) {
			Number tmp(x);
			tmp += i;
			retval *= tmp;
		}
		return retval;
	}

	/// Rising factorial, double-precision argument.
	/**
	 * Will use piranha::generic_r_factorial internally.
	 * @throws value_error if i is negative.
	 */
	inline double r_factorial(const double &x, const int &i)
	{
		return generic_r_factorial(x,i);
	}

	/// Rising factorial, complex double-precision argument.
	/**
	 * Will use piranha::generic_r_factorial internally.
	 * @throws value_error if i is negative.
	 */
	inline std::complex<double> r_factorial(const std::complex<double> &c, const int &i)
	{
		return generic_r_factorial(c,i);
	}

	template <class Number>
	inline Number generic_f_factorial(const Number &x, const int &n)
	{
		factorial_check(n);
		Number retval(1);
		for (int i = 0; i < n; ++i) {
			Number tmp(x);
			tmp -= i;
			retval *= tmp;
		}
		return retval;
	}

	/// Falling factorial, double-precision argument.
	/**
	 * Will use piranha::generic_f_factorial internally.
	 * @throws value_error if i is negative.
	 */
	inline double f_factorial(const double &x, const int &i)
	{
		return generic_f_factorial(x,i);
	}

	/// Falling factorial, complex double-precision argument.
	/**
	 * Will use piranha::generic_f_factorial internally.
	 * @throws value_error if i is negative.
	 */
	inline std::complex<double> f_factorial(const std::complex<double> &c, const int &i)
	{
		return generic_f_factorial(c,i);
	}

	template <class Number>
	inline Number generic_choose(const Number &n, const int &k)
	{
		Number retval(0);
		if (k < 0) {
			return retval;
		}
		retval = 1;
		for (int i = 1; i <= k; ++i) {
			Number tmp(n);
			tmp -= k;
			tmp += i;
			retval *= tmp;
			retval /= i;
		}
		return retval;
	}

	/// Binomial coefficient.
	/**
	 * Will use internally Boost's implementation.
	 */
	inline double choose(const int &n_, const int &k)
	{
		const int multiplier = (n_ >= 0) ? 1 : cs_phase(k), n = (n_ >= 0) ? n_ : k - n_ - 1;
		if (k < 0 || k > n) {
			return 0.;
		}
		piranha_assert(n >= 0);
		return boost::math::binomial_coefficient<double>(
			static_cast<unsigned>(n),static_cast<unsigned>(k)) * multiplier;
	}

	/// Binomial coefficient.
	/**
	 * Will use piranha::generic_choose internally.
	 */
	inline double choose(const double &x, const int &k)
	{
		return generic_choose(x,k);
	}

	/// Binomial coefficient.
	/**
	 * Will use piranha::generic_choose internally.
	 */
	inline std::complex<double> choose(const std::complex<double> &c, const int &k)
	{
		return generic_choose(c,k);
	}

	/// Binomial coefficient.
	/**
	 * Will use piranha::generic_choose internally.
	 */
	inline std::complex<double> choose(const std::complex<int> &c, const int &k)
	{
		return generic_choose(std::complex<double>(c.real(),c.imag()),k);
	}

	/// Calculate complex exponential of n*pi/2.
	inline std::complex<double> einpi2(const int &n)
	{
		if (n & 1) {
			if ((n - 1) & 3) {
				return std::complex<double>(0,-1);
			} else {
				return std::complex<double>(0,1);
			}
		} else {
			if (n & 3) {
				return std::complex<double>(-1,0);
			} else {
				return std::complex<double>(1,0);
			}
		}
	}

	/// Bessel function of the first kind, integer order, real double-precision argument.
	inline double besselJ(const int &order, const double &arg)
	{
		return boost::math::cyl_bessel_j(order, arg);
	}

	/// First derivative of the Bessel function of the first kind, integer order, double-precision argument.
	inline double dbesselJ(const int &order, const double &arg)
	{
		return (besselJ(order - 1,arg) - besselJ(order + 1, arg)) / 2.;
	}

	/// Bessel function of the first kind, integer order, complex double-precision argument.
	inline std::complex<double> besselJ(const int &order, const std::complex<double> &arg)
	{
		if (arg.imag() != 0) {
			piranha_throw(value_error,"only complex values with null imaginary part can be used as arguments "
				"for besselJ");
		}
		return std::complex<double>(boost::math::cyl_bessel_j(order, arg.real()));
	}

	/// First derivative ot the Bessel function of the first kind, integer order, complex double-precision argument.
	inline std::complex<double> dbesselJ(const int &order, const std::complex<double> &arg)
	{
		if (arg.imag() != 0) {
			piranha_throw(value_error,"only complex values with null imaginary part can be used as arguments "
				"for dbesselJ");
		}
		return std::complex<double>(dbesselJ(order, arg.real()));
	}

	/// Modified Bessel function of the first kind, integer order.
	inline double besselI(const int &order, const double &arg)
	{
		return boost::math::cyl_bessel_i(order, arg);
	}

	/// Associated Legendre function Pnm.
	inline double legendrePnm(const int &n, const int &m, const double &arg)
	{
		return boost::math::legendre_p(n, m, arg);
	}

	/// Legendre polynomial Pn.
	inline double legendrePn(const int &n, const double &arg)
	{
		return boost::math::legendre_p(n, arg);
	}

	inline std::complex<double> Ynm(const int &n, const int &m, const double &theta, const std::complex<double> &ei_phi,
		const std::complex<double> &)
	{
		std::complex<double> retval(std::pow(ei_phi,m));
		retval *= legendrePnm(n,m,std::cos(theta));
		return retval;
	}

	/// Non-normalised spherical harmonic.
	inline std::complex<double> Ynm(const int &n, const int &m, const double &theta, const double &phi)
	{
		return Ynm(n,m,theta,std::polar(1.,phi),std::complex<double>());
	}

	/// Non-normalised rotated spherical harmonic.
	inline std::complex<double> Ynm(const int &n_, const int &m_, const double &theta, const std::complex<double> &ei_phi,
		const std::complex<double> &, const double &alpha, const double &beta, const double &gamma)
	{
		// Let's fix negative n and/or m.
		int n(n_), m(std::abs(m_));
		std::complex<double> retval(1.,0.);
		if (n_ < 0) {
			n = -n_-1;
		}
		if (n == 0 && m == 0) {
			return retval;
		}
		retval = std::complex<double>(0.,0.);
		if (m > n) {
			return retval;
		}
		piranha_assert(n >= m && n >= 0 && m >= 0);
		std::complex<double> factor(std::polar(1.,-gamma*m));
		factor *= einpi2(-m);
		factor *= factorial(n+m);
		for (int k = -n; k <=n; ++k) {
			std::complex<double> tmp = std::pow(ei_phi,k) * std::polar(1.,-alpha*k);
			tmp *= einpi2(k);
			tmp *= legendrePnm(n,k,std::cos(theta));
			tmp *= factorial(n-k);
			double tmp2 = 0;
			for (int t = std::max<int>(0,k-m); t <= std::min<int>(n-m,n+k); ++t) {
				double tmp3 = std::pow(std::sin(beta/2.),m-k+2*t) * std::pow(std::cos(beta/2.),2*n-m+k-2*t);
				tmp3 *= cs_phase(t);
				tmp3 /= factorial(t)*factorial(n+k-t)*factorial(n-m-t)*factorial(m-k+t);
				tmp2 += tmp3;
			}
			tmp *= tmp2;
			retval += tmp;
		}
		return retval * factor;
	}

	inline std::complex<double> Ynm(const int &n, const int &m, const double &theta, const double &phi,
		const double &alpha, const double &beta, const double &gamma)
	{
		return Ynm(n,m,theta,std::polar(1.,phi),std::complex<double>(),alpha,beta,gamma);
	}

	// Root functions.

	/// Root function.
	/**
	 * @throws value_error if n <= 0 or if x < 0 and n is not unitary.
	 */
	inline double root(const int &n, const double &x)
	{
		if (n <= 0) {
			piranha_throw(value_error,"root order must be strictly positive");
		}
		if (x < 0 && n > 1) {
			piranha_throw(value_error,"cannot calculate the root of a negative real number");
		}
		return std::pow(x, 1. / (double)n);
	}

	/// Root function, complex version.
	/**
	 * @throws value_error if n <= 0.
	 */
	inline std::complex<double> root(const int &n, const std::complex<double> &c)
	{
		if (n <= 0) {
			piranha_throw(value_error,"root order must be strictly positive");
		}
		return std::pow(c, 1. / (double)n);
	}

	/// Multiply-accumulate.
	/**
	 * Equivalent to x += y * z.
	 */
	template <class T, class U, class V>
	inline void multiply_accumulate(T &x, const U &y, const V &z)
	{
		x += y * z;
	}

	/// Check if input argument can be represented exactly as an int type.
	/**
	 * Returns true.
	 */
	inline bool is_integer(const int &)
	{
		return true;
	}

	/// Check if input argument can be represented exactly as an int type.
	__PIRANHA_VISIBLE bool is_integer(const double &);
}

#endif
