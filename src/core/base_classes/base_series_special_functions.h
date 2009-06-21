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
#include <vector>

#include "../exceptions.h"
#include "../math.h"
#include "../mp.h"
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
			template <class ArgsTuple>
			Derived base_hyperF(const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list, const int &iter_limit, const ArgsTuple &args_tuple) const
			{
				Derived retval;
				// If one of the elements in b_list is a non-positive integer, a division by zero will occur.
				hyperF_b_check(b_list);
				// Cache values.
				const size_t a_size = a_list.size(), b_size = b_list.size();
				// HyperF of a null series will always be equal to 1.
				if (derived_const_cast->empty()) {
					retval.base_add(1,args_tuple);
					return retval;
				}
				// If one of the elements in a_list is zero, hyperF will always be 1.
				for (size_t i = 0; i < a_size; ++i) {
					if (a_list[i] == 0) {
						retval.base_add(1,args_tuple);
						return retval;
					}
				}
				// Establish series limit, using the following strategy...
				size_t iter_;
				try {
					// First we try to respect the given iter_limit, if any,
					// otherwise we call the truncator to see what it has to say.
					iter_ = (iter_limit >= 0) ? iter_limit : derived_const_cast->psi_(0,1,args_tuple);
				} catch (const value_error &) {
					// If the truncator fails to establish a limit, we go to see into a_list and b_list
					// whether there are negative integer numbers.
					mp_rational const *min_int = 0;
					for (size_t i = 0; i < a_size; ++i) {
						if (a_list[i] < 0 && a_list[i].get_den() == 1) {
							if (!min_int || a_list[i] < *min_int) {
								min_int = &a_list[i];
							}
						}
					}
					if (!min_int) {
						piranha_throw(value_error,"could not establish a value for the series limit in hyperF: "
							"no limit was provided, the truncator failed to establish a limit and "
							"no negative integer numbers were found in a_list");
					}
					piranha_assert(*min_int < 0);
					iter_ = -(((*min_int) - 1).to_int());
				}
				const size_t iter = iter_;
				// If no iterations are needed, return an empty series.
				if (!iter) {
					return retval;
				}
				retval.base_add(1,args_tuple);
				Derived tmp;
				tmp.base_add(1,args_tuple);
				for (size_t i = 1; i < iter; ++i) {
					for (size_t j = 0; j < a_size; ++j) {
						tmp.base_mult_by(a_list[j] + (int)(i - 1),args_tuple);
					}
					for (size_t j = 0; j < b_size; ++j) {
						tmp.base_divide_by(b_list[j] + (int)(i - 1),args_tuple);
					}
					tmp.base_divide_by(i,args_tuple);
					tmp.base_mult_by(*derived_const_cast,args_tuple);
					retval.base_add(tmp,args_tuple);
				}
				return retval;
			}
			template <class ArgsTuple>
			Derived base_dhyperF(const int &n, const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list, const int &iter_limit, const ArgsTuple &args_tuple) const
			{
				if (n < 0) {
					piranha_throw(value_error,"please enter a non-negative order of differentiation");
				}
				// If one of the elements in b_list is a non-positive integer, a division by zero will occur.
				hyperF_b_check(b_list);
				// Compute the outside factor.
				mp_rational factor(1);
				const size_t a_size = a_list.size(), b_size = b_list.size();
				for (size_t i = 0; i < a_size; ++i) {
					factor *= a_list[i].r_factorial(n);
				}
				for (size_t i = 0; i < b_size; ++i) {
					factor /= b_list[i].r_factorial(n);
				}
				// Create the new a_list and b_list.
				std::vector<mp_rational> new_a(a_list), new_b(b_list);
				for (size_t i = 0; i < a_size; ++i) {
					new_a[i] += n;
				}
				for (size_t i = 0; i < b_size; ++i) {
					new_b[i] += n;
				}
				Derived retval(base_hyperF(new_a,new_b,iter_limit,args_tuple));
				retval.base_mult_by(factor,args_tuple);
				return retval;
			}
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
				// Note: copy ctor of named derived class is not a problem, since it does not trim.
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
		private:
			static void hyperF_b_check(const std::vector<mp_rational> &b_list)
			{
				const size_t b_size = b_list.size();
				for (size_t i = 0; i < b_size; ++i) {
					if (b_list[i] <= 0 && b_list[i].get_den() == 1) {
						piranha_throw(zero_division_error,"b_list in hyperF contains a non-positive integer");
					}
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
