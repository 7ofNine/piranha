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

#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../math.h"
#include "../mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class base_series_special_functions
	{
		//protected:
		public:

			template <class ArgsTuple>
			Derived base_hyperF(const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list, const int &iter_limit, const ArgsTuple &argsTuple) const
			{
				Derived retval;
				// If one of the elements in b_list is a non-positive integer, a division by zero will occur.
				hyperF_b_check(b_list);
				// Cache values.
				const std::size_t a_size = a_list.size(), b_size = b_list.size();
				// HyperF of a null series will always be equal to 1.
				if (derived_const_cast->empty()) 
				{
					retval.base_add(1,argsTuple);
					return retval;
				}
				// If one of the elements in a_list is zero, hyperF will always be 1.
				for (std::size_t i = 0; i < a_size; ++i) 
				{
					if (a_list[i] == 0) 
					{
						retval.base_add(1,argsTuple);
						return retval;
					}
				}
				// Establish series limit, using the following strategy...
				std::size_t iter_;
				try {
					// First we try to respect the given iter_limit, if any,
					// otherwise we call the truncator to see what it has to say.
					iter_ = (iter_limit >= 0) ? iter_limit : derived_const_cast->psi_(0,1,argsTuple);
				} catch (const value_error &) 
				{
					// If the truncator fails to establish a limit, we go to see into a_list and b_list
					// whether there are negative integer numbers.
					mp_rational const *min_int = 0;
					for (std::size_t i = 0; i < a_size; ++i) 
					{
						if (a_list[i] < 0 && a_list[i].get_den() == 1) 
						{
							if (!min_int || a_list[i] < *min_int) 
							{
								min_int = &a_list[i];
							}
						}
					}

					if (!min_int) 
					{
						PIRANHA_THROW(value_error,"could not establish a value for the series limit in hyperF: "
							"no explicit limit was provided, the truncator failed to establish a limit and "
							"no negative integer numbers were found in a_list");
					}
					PIRANHA_ASSERT(*min_int < 0);
					iter_ = -(((*min_int) - 1).to_int());
				}

				const std::size_t iter = iter_;
				// If no iterations are needed, return an empty series.
				if (!iter) 
				{
					return retval;
				}

				retval.base_add(1,argsTuple);
				Derived tmp;
				tmp.base_add(1,argsTuple);
				for (std::size_t i = 1; i < iter; ++i) 
				{
					for (std::size_t j = 0; j < a_size; ++j) 
					{
						tmp.base_mult_by(a_list[j] + (int)(i - 1), argsTuple);
					}

					for (std::size_t j = 0; j < b_size; ++j) 
					{
						tmp.base_divide_by(b_list[j] + (int)(i - 1), argsTuple);
					}

					tmp.base_divide_by(boost::numeric_cast<int>(i), argsTuple);
					tmp.base_mult_by(*derived_const_cast,argsTuple);
					retval.base_add(tmp,argsTuple);
				}

				return retval;
			}


			template <class ArgsTuple>
			Derived base_dhyperF(const int &n, const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list, const int &iter_limit, const ArgsTuple &argsTuple) const
			{
				if (n < 0) 
				{
					PIRANHA_THROW(value_error,"please enter a non-negative order of differentiation");
				}

				// If one of the elements in b_list is a non-positive integer, a division by zero will occur.
				hyperF_b_check(b_list);
				// Compute the outside factor.
				mp_rational factor(1);
				const std::size_t a_size = a_list.size(), b_size = b_list.size();
				for (std::size_t i = 0; i < a_size; ++i) 
				{
					factor *= a_list[i].r_factorial(n);
				}

				for (std::size_t i = 0; i < b_size; ++i) 
				{
					factor /= b_list[i].r_factorial(n);
				}

				// Create the new a_list and b_list.
				std::vector<mp_rational> new_a(a_list), new_b(b_list);
				for (std::size_t i = 0; i < a_size; ++i) 
				{
					new_a[i] += n;
				}

				for (std::size_t i = 0; i < b_size; ++i) {
					new_b[i] += n;
				}

				Derived retval(base_hyperF(new_a, new_b, iter_limit, argsTuple));
				retval.base_mult_by(factor, argsTuple);
				return retval;
			}


			/// Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived base_besselJ(const int &order_, const ArgsTuple &argsTuple) const 
			{
				Derived retval;
				// Dispatch besselJ to coefficient if series consists of a single coefficient.
				if (derived_const_cast->isSingleCf()) 
				{
					typedef typename Derived::TermType term_type;
					retval.insert(term_type(derived_const_cast->begin()->cf.besselJ(order_, argsTuple),
						typename term_type::key_type()), argsTuple);
					return retval;
				}

				// Take care of negative order.
				const int order = (order_ >= 0) ? order_ : -order_;
				// Calculate this/2.
				Derived x_2(*derived_const_cast);
				x_2.base_divide_by(2, argsTuple);
				retval = Derived(x_2).base_mult_by(-1, argsTuple).base_mult_by(x_2, argsTuple)
					.base_hyperF(std::vector<mp_rational>(), std::vector<mp_rational>((std::size_t)1, mp_rational(order) + 1), -1, argsTuple);
				retval.base_mult_by(x_2.base_pow(order, argsTuple), argsTuple);
				retval.base_divide_by(mp_integer(order).factorial(), argsTuple);
				if (order_ < 0) 
				{
					retval.base_mult_by(cs_phase(order_), argsTuple);
				}
				return retval;
			}


			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived base_dbesselJ(const int &order_, const ArgsTuple &argsTuple) const
			{
				// NOTE: we don't user hyperF here because direct series expansion is faster.
				const int order = (order_ >= 0) ? order_ : -order_;
				Derived retval;
				if (!order) 
				{
					// Exploit the fact that for order == 0, derivative is -J_1.
					retval = base_besselJ(1, argsTuple);
					retval.base_mult_by(-1, argsTuple);
					return retval;
				}

				// Get the expansion limit from the truncator.
				std::size_t limit_;
				try {
					limit_ = derived_const_cast->psi_(order - 1, 2, argsTuple);
				} catch (const value_error &ve) 
				{
					PIRANHA_THROW(value_error,std::string("series is unsuitable as argument of the derivative of "
						"Bessel function of the first kind.\nThe reported error is: ")
						+ ve.what());
				}

				const std::size_t limit = limit_;
				if (!limit) 
				{
					return retval;
				}
				// Now we buid the starting point of the power series expansion.
				retval = *derived_const_cast;
				retval.base_divide_by(2, argsTuple);
				// This will be used later.
				Derived square_x2(retval);
				square_x2.base_mult_by(square_x2, argsTuple);
				retval = retval.base_pow(order - 1, argsTuple);
				retval.base_divide_by(mp_integer(order).factorial(), argsTuple);
				retval.base_mult_by(order, argsTuple);
				retval.base_divide_by(2, argsTuple);
				// Now let's proceed to the bulk of the power series expansion.
				Derived tmp(retval);
				for (int i = 1; i < (int)limit; ++i) 
				{
					tmp.base_mult_by((-1) * (mp_integer(order) + 2 * mp_integer(i)), argsTuple);
					tmp.base_divide_by((mp_integer(i) * (mp_integer(i) + order)) * (mp_integer(order) + 2 * (mp_integer(i) - 1)), argsTuple);
					tmp.base_mult_by(square_x2, argsTuple);
					retval.base_add(tmp, argsTuple);
				}
				return retval;
			}


			/// Bessel function of the first kind of integer order divided by its argument**m.
			template <class ArgsTuple>
			Derived base_besselJ_div_m(const int &order_, const int &m, const ArgsTuple &argsTuple) const
			{
				// Take care of negative order.
				const int order = (order_ >= 0) ? order_ : -order_;
				if (order < m) 
				{
					PIRANHA_THROW(zero_division_error,"absolute value of order must not be smaller than m");
				}

				Derived retval;
				// Calculate this/2.
				Derived x_2(*derived_const_cast);
				x_2.base_divide_by(2, argsTuple);

				retval = Derived(x_2).base_mult_by(-1, argsTuple).base_mult_by(x_2, argsTuple)
					.base_hyperF(std::vector<mp_rational>(), std::vector<mp_rational>((std::size_t)1, mp_rational(order) + 1), -1, argsTuple);
				retval.base_mult_by(derived_const_cast->base_pow(order - m, argsTuple), argsTuple);
				retval.base_divide_by(mp_integer(order).factorial() * mp_integer(2).pow(order), argsTuple);

				if (order_ < 0) 
				{
					retval.base_mult_by(cs_phase(order_), argsTuple);
				}

				return retval;
			}

		private:

			static void hyperF_b_check(const std::vector<mp_rational> &b_list)
			{
				const std::size_t b_size = b_list.size();
				for (std::size_t i = 0; i < b_size; ++i) 
				{
					if (b_list[i] <= 0 && b_list[i].get_den() == 1) 
					{
						PIRANHA_THROW(zero_division_error, "b_list in hyperF contains a non-positive integer");
					}
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
