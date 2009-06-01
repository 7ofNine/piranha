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
#include "../math.h"
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	struct base_series_special_functions {};

	template <class Derived>
	class toolbox<base_series_special_functions<Derived> >
	{
		protected:
			/// Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived base_besselJ(const int &order_, const ArgsTuple &args_tuple) const {
				Derived retval;
				// Special case of empty series. It must be handled here before the truncator is called.
				if (derived_const_cast->empty()) {
					if (!order_) {
						retval.base_add(1,args_tuple);
					}
					return retval;
				}
				// Dispatch besselJ to coefficient if series consists of a single coefficient.
				if (derived_const_cast->is_single_cf()) {
					typedef typename Derived::term_type term_type;
					retval.insert(term_type(derived_const_cast->begin()->m_cf.besselJ(order_,args_tuple),
						typename term_type::key_type()),args_tuple);
					return retval;
				}
				const int order = (order_ >= 0) ? order_ : -order_;
				// Get the expansion limit from the truncator.
				size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order, 2, args_tuple);
				} catch (const value_error &ve) {
					piranha_throw(value_error,std::string("series is unsuitable as argument of "
						"Bessel function of the first kind.\nThe reported error is: ") + ve.what());
				}
				const size_t limit = limit_;
				if (limit == 0) {
					return retval;
				}
				// Now we buid the starting point of the power series expansion of Jn.
				retval = *derived_const_cast;
				retval.base_divide_by(2, args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.base_mult_by(square_x2, args_tuple);
				retval = retval.base_pow(order, args_tuple);
				retval.base_divide_by(factorial(order),args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn.
				Derived tmp(retval);
				for (size_t i = 1; i < limit; ++i) {
					tmp.base_mult_by(-1, args_tuple);
					tmp.base_divide_by(i * ((double)i + order), args_tuple);
					tmp.base_mult_by(square_x2, args_tuple);
					retval.base_add(tmp, args_tuple);
				}
				if (order_ < 0) {
					retval.base_mult_by(cs_phase(order_), args_tuple);
				}
				return retval;
			}
			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			// TODO: polish and make it less restrictive. Use factorial().
			template <class ArgsTuple>
			Derived base_dbesselJ(const int &order, const ArgsTuple &args_tuple) const {
				if (order < 1) {
					piranha_throw(value_error,"partial derivative of Bessel function of the first kind is implemented only for "
									 "strictly positive orders");
				}
				// Get the expansion limit from the truncator.
				size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order - 1, 2, args_tuple);
				} catch (const value_error &ve) {
					piranha_throw(value_error,std::string("series is unsuitable as argument of the derivative of "
						"Bessel function of the first kind.\nThe reported error is: ")
						+ ve.what());
				}
				const size_t limit = limit_;
				Derived retval;
				if (limit == 0) {
					return retval;
				}
				// Now we buid the starting point of the power series expansion of Jn.
				retval = *derived_const_cast;
				retval.base_divide_by(2, args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.base_mult_by(square_x2, args_tuple);
				retval = retval.base_pow(order - 1, args_tuple);
				for (size_t i = 0; i < (size_t)order; ++i) {
					retval.base_divide_by(i + 1, args_tuple);
				}
				retval.base_mult_by(order, args_tuple);
				retval.base_divide_by(2, args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn.
				Derived tmp(retval);
				for (size_t i = 1; i < limit; ++i) {
					tmp.base_mult_by((-1) * (order + 2 * (double)i), args_tuple);
					tmp.base_divide_by(((double)i * ((double)i + order)) * (order + 2 * ((double)i - 1)), args_tuple);
					tmp.base_mult_by(square_x2, args_tuple);
					retval.base_add(tmp, args_tuple);
				}
				return retval;
			}
			/// Bessel function of the first kind of integer order divided by its argument**m.
			template <class ArgsTuple>
			Derived base_besselJ_div_m(const int &order_, const int &m,
				const ArgsTuple &args_tuple) const
			{
				Derived retval;
				// Let's take care of negative order.
				const int order = (order_ >= 0) ? order_ : -order_;
				// Special case of empty series.
				if (derived_const_cast->empty()) {
					if (order < m) {
						piranha_throw(zero_division_error,"besselJ_div_m with order smaller than m has a pole in the origin");
					} else if (order == m) {
						retval = retval.base_add(2,args_tuple).base_pow(-order,args_tuple);
						retval.base_divide_by(factorial(order),args_tuple);
					}
					return retval;
				}
				// Get the expansion limit from the truncator.
				size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order - m, 2, args_tuple);
				} catch (const value_error &ve) {
					piranha_throw(value_error,std::string("series is unsuitable as argument of Bessel function "
						"of the first kind divided by its argument raised to unsigned integer power.\n"
						"The reported error is: ") + ve.what());
				}
				const size_t limit = limit_;
				if (limit == 0) {
					return retval;
				}
				// Now we build the starting point of the power series expansion of Jn/x**m.
				retval = *derived_const_cast;
				retval.base_divide_by(2, args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.base_mult_by(square_x2, args_tuple);
				retval = retval.base_pow(order - m, args_tuple);
				retval.base_divide_by(factorial(order),args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn/x**m.
				Derived tmp(retval);
				for (size_t i = 1; i < limit; ++i) {
					tmp.base_mult_by(-1, args_tuple);
					tmp.base_divide_by(i * ((double)i + order), args_tuple);
					tmp.base_mult_by(square_x2, args_tuple);
					retval.base_add(tmp, args_tuple);
				}
				// m2 will be the external divisor by 2**m.
				Derived m2;
				m2 = m2.base_add(2,args_tuple).base_pow(-m,args_tuple);
				retval.base_mult_by(m2, args_tuple);
				if (order_ < 0) {
					retval.base_mult_by(cs_phase(order_), args_tuple);
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
