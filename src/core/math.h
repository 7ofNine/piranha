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
#include "integer_typedefs.h"
#include "p_assert.h"

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

	template <class T>
	inline max_fast_int sign(const T &x)
	{
		if (x >= 0) {
			return static_cast<max_fast_int>(1);
		} else {
			return static_cast<max_fast_int>(-1);
		}
	}

	/// Condon-Shortley phase.
	inline max_fast_int cs_phase(const max_fast_int &n)
	{
		if (n & 1) {
			return static_cast<max_fast_int>(-1);
		} else {
			return static_cast<max_fast_int>(1);
		}
	}

	/// Swap two integers without using extra storage.
	/**
	 * The two integers must not be in the same memory location.
	 */
	template <class Integer>
	inline void int_swap(Integer &a, Integer &b)
	{
		p_assert(&a != &b);
		(((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)));
	}

	template <class T>
	static inline void factorial_check(const T &x) {
		if (x < 0) {
			throw unsuitable("Please use a non-negative integer as argument for factorials.");
		}
	}

	/// Factorial.
	inline double factorial(const int &i) {
		factorial_check(i);
		return boost::math::factorial<double>(static_cast<unsigned>(i));
	}

	/// Double factorial.
	inline double double_factorial(const int &i) {
		factorial_check(i);
		return boost::math::double_factorial<double>(static_cast<unsigned>(i));
	}

	/// Binomial coefficient.
	inline double choose(const max_fast_int &n, const max_fast_int &k)
	{
		if (n < 0 || k < 0 || k > n) {
			throw unsuitable("Invalid input values for binomial coefficient.");
		}
		return boost::math::binomial_coefficient<double>(
			static_cast<unsigned>(n),static_cast<unsigned>(k));
	}

	/// Calculate complex exponential of n*pi/2.
	inline std::complex<max_fast_int> einpi2(const max_fast_int &n)
	{
		if (n & 1) {
			if ((n-1) & 3) {
				return std::complex<max_fast_int>(static_cast<max_fast_int>(0),static_cast<max_fast_int>(-1));
			} else {
				return std::complex<max_fast_int>(static_cast<max_fast_int>(0),static_cast<max_fast_int>(1));
			}
		} else {
			if (n & 3) {
				return std::complex<max_fast_int>(static_cast<max_fast_int>(-1),static_cast<max_fast_int>(0));
			} else {
				return std::complex<max_fast_int>(static_cast<max_fast_int>(1),static_cast<max_fast_int>(0));
			}
		}
	}

	/// Bessel function of the first kind, integer order.
	inline double besselJ(const max_fast_int &order, const double &arg)
	{
		return boost::math::cyl_bessel_j(order, arg);
	}

	/// Modified Bessel function of the first kind, integer order.
	inline double besselI(const max_fast_int &order, const double &arg)
	{
		return boost::math::cyl_bessel_i(order, arg);
	}

	/// Associated Legendre function Pnm.
	inline double Pnm(const int &n, const int &m, const double &arg)
	{
		return boost::math::legendre_p(n, m, arg);
	}

	/// Legendre polynomial Pn.
	inline double Pn(const int &n, const double &arg)
	{
		return boost::math::legendre_p(n, arg);
	}

	inline std::complex<double> Ynm(const int &n, const int &m, const double &theta, const std::complex<double> &ei_phi,
		const std::complex<double> &)
	{
		std::complex<double> retval(std::pow(ei_phi,m));
		retval *= Pnm(n,m,std::cos(theta));
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
		p_assert(n >= m && n >= 0 && m >= 0);
		std::complex<double> factor(std::polar(1.,-gamma*m));
		factor *= einpi2(-m);
		factor *= factorial(n+m);
		for (int k = -n; k <=n; ++k) {
			std::complex<double> tmp = std::pow(ei_phi,k) * std::polar(1.,-alpha*k);
			tmp *= einpi2(k);
			tmp *= Pnm(n,k,std::cos(theta));
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
}

#endif
