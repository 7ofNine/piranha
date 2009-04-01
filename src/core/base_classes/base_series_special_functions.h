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
			Derived besselJ_(const max_fast_int &order_, const ArgsTuple &args_tuple) const {
				Derived retval;
				// Special case of empty series. It must be handled here before the truncator is called.
				if (derived_const_cast->empty()) {
					if (!order_) {
						retval = Derived(max_fast_int(1),args_tuple);
					}
					return retval;
				}
				// Dispatch besselJ to coefficient if series consists of a single coefficient.
				if (derived_const_cast->is_single_cf()) {
					typedef typename Derived::term_type term_type;
					retval.insert(term_type(derived_const_cast->begin()->m_cf.besselJ_(order_,args_tuple),
						typename term_type::key_type()),args_tuple);
					return retval;
				}
				const max_fast_int order = (order_ >= 0) ? order_ : -order_;
				// Get the expansion limit from the truncator.
				size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order, 2, args_tuple);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable as argument of "
						"Bessel function of the first kind.\nThe reported error is: ") + u.what());
				}
				const size_t limit = limit_;
				// Now we buid the starting point of the power series expansion of Jn.
				retval = *derived_const_cast;
				retval.divide_by(max_fast_int(2), args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.mult_by(square_x2, args_tuple);
				retval = retval.pow_(order, args_tuple);
				retval.mult_by(Derived::factorial_(order,args_tuple).inv_(args_tuple),args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn.
				Derived tmp(retval);
				for (size_t i = 1; i < limit; ++i) {
					tmp.mult_by(max_fast_int(-1), args_tuple);
					tmp.divide_by(max_fast_int(i*(i + order)), args_tuple);
					tmp.mult_by(square_x2, args_tuple);
					retval.add(tmp, args_tuple);
				}
				if (order_ < 0) {
					retval.mult_by(cs_phase(order_), args_tuple);
				}
				return retval;
			}
			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			// TODO: polish and make it less restrictive. Use factorial().
			template <class ArgsTuple>
			Derived dbesselJ_(const max_fast_int &order, const ArgsTuple &args_tuple) const {
				if (order < 1) {
					throw unsuitable("Partial derivative of Bessel function of the first kind is implemented only for "
									 "strictly positive orders.");
				}
				// Get the expansion limit from the truncator.
				size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order - 1, 2, args_tuple);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable as argument of the derivative of "
												 "Bessel function of the first kind.\nThe reported error is: ")
									 + u.what());
				}
				const size_t limit = limit_;
				// Now we buid the starting point of the power series expansion of Jn.
				Derived retval(*derived_const_cast);
				retval.divide_by(max_fast_int(2), args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.mult_by(square_x2, args_tuple);
				retval = retval.pow_(order - 1, args_tuple);
				for (size_t i = 0; i < (size_t)order; ++i) {
					retval.divide_by((max_fast_int)(i + 1), args_tuple);
				}
				retval.mult_by(order, args_tuple);
				retval.divide_by((max_fast_int)2, args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn.
				Derived tmp(retval);
				for (size_t i = 1; i < limit; ++i) {
					tmp.mult_by((max_fast_int)(-1)*(max_fast_int)(order + 2*i), args_tuple);
					tmp.divide_by((max_fast_int)(i*(i + order))*(max_fast_int)(order + 2*(i - 1)), args_tuple);
					tmp.mult_by(square_x2, args_tuple);
					retval.add(tmp, args_tuple);
				}
				return retval;
			}
			/// Bessel function of the first kind of integer order divided by its argument**m.
			template <class ArgsTuple>
			Derived besselJ_div_m_(const max_fast_int &order_, const max_fast_int &m,
				const ArgsTuple &args_tuple) const {
				Derived retval;
				// Let's take care of negative order.
				const max_fast_int order = (order_ >= 0) ? order_ : -order_;
				// Special case of empty series.
				if (derived_const_cast->empty()) {
					if (order < m) {
						throw division_by_zero();
					} else if (order == m) {
						retval = Derived(max_fast_int(2),args_tuple).pow_(-order,args_tuple);
						retval.mult_by(Derived::factorial_(order,args_tuple).pow_(-1,args_tuple),args_tuple);
					}
					return retval;
				}
				// Get the expansion limit from the truncator.
				size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order - m, 2, args_tuple);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable as argument of Bessel function "
						"of the first kind divided by its argument raised to unsigned integer power.\n"
						"The reported error is: ") + u.what());
				}
				const size_t limit = limit_;
				// Now we build the starting point of the power series expansion of Jn/x**m.
				retval = *derived_const_cast;
				retval.divide_by(max_fast_int(2), args_tuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.mult_by(square_x2, args_tuple);
				retval = retval.pow_(order - m, args_tuple);
				retval.mult_by(Derived::factorial_(order,args_tuple).inv_(args_tuple),args_tuple);
				// Now let's proceed to the bulk of the power series expansion for Jn/x**m.
				Derived tmp(retval);
				for (size_t i = 1; i < limit; ++i) {
					tmp.mult_by(max_fast_int(-1), args_tuple);
					tmp.divide_by(max_fast_int(i*(i + order)), args_tuple);
					tmp.mult_by(square_x2, args_tuple);
					retval.add(tmp, args_tuple);
				}
				// m2 will be the external divisor by 2**m.
				Derived m2 = Derived(max_fast_int(2),args_tuple).pow_(-m,args_tuple);
				retval.mult_by(m2, args_tuple);
				if (order_ < 0) {
					retval.mult_by(cs_phase(order_), args_tuple);
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
