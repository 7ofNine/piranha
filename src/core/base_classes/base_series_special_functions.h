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

#ifndef PIRANHA_BASE_SERIES_SPECIAL_FUNCTIONS_H
#define PIRANHA_BASE_SERIES_SPECIAL_FUNCTIONS_H

#include <string>

#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class base_series_special_functions
	{
		public:
			/// Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived besselJ(const max_fast_int &order_, const ArgsTuple &args_tuple) const {
				size_t order;
				max_fast_int multiplier = 1;
				if (order_ >= 0) {
					order = order_;
				} else {
					order = static_cast<size_t>(-order_);
					if ((order & 1) != 0) {
						multiplier = -1;
					}
				}
				// Get the expansion limit from the truncator.
				size_t limit;
				try {
					limit = Derived::multiplier_type::truncator_type::power_series_limit(
								*derived_const_cast, args_tuple, order, 2);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable as argument of Bessel function of the first kind.\n"
												 "The reported error is: ")
									 + u.what());
				}
				// Now we buid the starting point of the power series expansion of Jn.
				Derived retval(*derived_const_cast);
				retval.divide_by((max_fast_int)2, args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.mult_by(square_x2, args_tuple);
				retval = retval.pow((max_fast_int)order, args_tuple);
				for (size_t i = 0; i < order; ++i) {
					retval.divide_by((max_fast_int)(i + 1), args_tuple);
				}
				// Now let's proceed to the bulk of the power series expansion for Jn.
				Derived tmp(retval);
				for (size_t i = 1; i <= limit; ++i) {
					tmp.mult_by((max_fast_int) - 1, args_tuple);
					tmp.divide_by((max_fast_int)(i*(i + order)), args_tuple);
					tmp.mult_by(square_x2, args_tuple);
					retval.add(tmp, args_tuple);
				}
				retval.mult_by(multiplier, args_tuple);
				return retval;
			}
			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived dbesselJ(const max_fast_int &order, const ArgsTuple &args_tuple) const {
				if (order < 1) {
					throw unsuitable("Partial derivative of Bessel function of the first kind is implemented only for "
									 "strictly positive orders.");
				}
				// Get the expansion limit from the truncator.
				size_t limit;
				try {
					limit = Derived::multiplier_type::truncator_type::power_series_limit(*derived_const_cast,
							args_tuple, order - 1, 2);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable as argument of the derivative of "
												 "Bessel function of the first kind.\nThe reported error is: ")
									 + u.what());
				}
				// Now we buid the starting point of the power series expansion of Jn.
				Derived retval(*derived_const_cast);
				retval.divide_by((max_fast_int)2, args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.mult_by(square_x2, args_tuple);
				retval = retval.pow(order - 1, args_tuple);
				for (size_t i = 0; i < (size_t)order; ++i) {
					retval.divide_by((max_fast_int)(i + 1), args_tuple);
				}
				retval.mult_by(order, args_tuple);
				retval.divide_by((max_fast_int)2, args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn.
				Derived tmp(retval);
				for (size_t i = 1; i <= limit; ++i) {
					tmp.mult_by((max_fast_int)(-1)*(max_fast_int)(order + 2*i), args_tuple);
					tmp.divide_by((max_fast_int)(i*(i + order))*(max_fast_int)(order + 2*(i - 1)), args_tuple);
					tmp.mult_by(square_x2, args_tuple);
					retval.add(tmp, args_tuple);
				}
				return retval;
			}
			/// Bessel function of the first kind of integer order divided by its argument.
			template <class ArgsTuple>
			Derived besselJ_div(const max_fast_int &order, const ArgsTuple &args_tuple) const {
				if (order < 1) {
					throw unsuitable("Bessel function of the first kind divided by its argument is implemented only for "
									 "strictly positive orders.");
				}
				// Get the expansion limit from the truncator.
				size_t limit;
				try {
					limit = Derived::multiplier_type::truncator_type::power_series_limit(*derived_const_cast,
							args_tuple, order - 1, 2);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable as argument of Bessel function of the first kind.\n"
												 "The reported error is: ")
									 + u.what());
				}
				// Now we buid the starting point of the power series expansion of Jn/x.
				Derived retval(*derived_const_cast);
				retval.divide_by((max_fast_int)2, args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.mult_by(square_x2, args_tuple);
				retval = retval.pow((max_fast_int)order - 1, args_tuple);
				for (size_t i = 0; i < (size_t)order; ++i) {
					retval.divide_by((max_fast_int)(i + 1), args_tuple);
				}
				// Now let's proceed to the bulk of the power series expansion for Jn/x.
				Derived tmp(retval);
				for (size_t i = 1; i <= limit; ++i) {
					tmp.mult_by((max_fast_int) - 1, args_tuple);
					tmp.divide_by((max_fast_int)(i*(i + order)), args_tuple);
					tmp.mult_by(square_x2, args_tuple);
					retval.add(tmp, args_tuple);
				}
				retval.divide_by((max_fast_int)2, args_tuple);
				return retval;
			}
			/// Legendre function of the first kind - Pnm(cos(theta)).
			/**
			* This implementation uses recurrence relations.
			*/
			template <class ArgsTuple>
			Derived Pnm(const max_fast_int &n_, const max_fast_int &m_, const ArgsTuple &args_tuple) const {
				// This is P00 right now.
				Derived retval(static_cast<max_fast_int>(1),args_tuple);
				max_fast_int n(n_), m(std::abs(m_));
				if (n_ < 0) {
					n = -n_-1;
				}
				if (n == 0 && m == 0) {
					return retval;
				}
				if (m > n) {
					retval = Derived();
					return retval;
				}
				p_assert(n >= m && n >= 0 && m >= 0);
				Derived P00(retval), old_Pnm, tmp1;
				const std::complex<Derived> expi_theta(derived_const_cast->complexp(args_tuple));
				const Derived cos_theta(expi_theta.real(args_tuple)), sin_theta(expi_theta.imag(args_tuple));
				max_fast_int i = 0;
				// Recursion to get from P_00 to P_mm.
				for (; i < m; ++i) {
					retval.mult_by(i*2-1,args_tuple);
					retval.mult_by(sin_theta,args_tuple);
				}
				p_assert(i == m);
				// Recursion to get from P_mm to P_nm (n>m).
				for (; i < n; ++i) {
					old_Pnm.mult_by(-m-i,args_tuple);
					old_Pnm.divide_by(-m+i+1,args_tuple);
					tmp1 = old_Pnm;
					old_Pnm = retval;
					retval.mult_by(i*2+1,args_tuple);
					retval.divide_by(-m+i+1,args_tuple);
					retval.mult_by(cos_theta,args_tuple);
					retval.add(tmp1,args_tuple);
				}
				if (m_ < 0) {
					for (max_fast_int j = n + m; j >= n - m + 1; --j) {
						retval.divide_by(j,args_tuple);
					}
					retval.mult_by(cs_phase(m),args_tuple);
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
