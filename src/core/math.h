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

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/static_assert.hpp>
#include <cmath>
#include <complex>

#include "exceptions.h"
#include "integer_typedefs.h"

namespace piranha
{
	/// Meta-programmed functor for the calculation of base-2 logarithm.
	/**
	 * Result is retrieved through the lg::value const member function.
	 */
	template <int N>
	struct lg {
		BOOST_STATIC_ASSERT(N > 0 && (N % 2) == 0);
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

	/// Non-normalised spherical harmonic.
	inline std::complex<double> Ynm(const int &n, const int &m, const double &theta, const double &phi)
	{
		std::complex<double> retval(std::polar(1.,phi*m));
		retval *= Pnm(n,m,std::cos(theta));
		return retval;
	}

	template <class T>
	static inline void factorial_check(const T &x) {
		if (x < 0) {
			throw unsuitable("Please use a non-negative integer as argument for factorials.");
		}
	}

	inline double factorial(const int &i) {
		factorial_check(i);
		return boost::math::factorial<double>(static_cast<unsigned>(i));
	}

	inline double double_factorial(const int &i) {
		factorial_check(i);
		return boost::math::double_factorial<double>(static_cast<unsigned>(i));
	}
}

#endif
