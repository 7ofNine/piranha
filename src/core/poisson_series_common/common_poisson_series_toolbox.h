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

#ifndef PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H

#include <cmath> // For nearbyint.

#include "../base_classes/series_math.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class common_poisson_series_toolbox
	{
		public:
			Derived sin() const {
				Derived retval;
				// In this case the Poisson series is logically equivalent to zero.
				if (derived_const_cast->empty()) {
					// Return an empty series, hence do nothing.
					;
				} else {
					linear_int_helper<false>(retval);
				}
				return retval;
			}
			Derived cos() const {
				typedef typename Derived::term_type term_type;
				typedef typename term_type::cf_type cf_type;
				typedef typename term_type::key_type key_type;
				Derived retval;
				// In this case the Poisson series is logically equivalent to zero.
				if (derived_const_cast->empty()) {
					// Return series logically equivalent to 1.
					retval.insert(term_type(cf_type((max_fast_int)1, retval.m_arguments), key_type()),
								  retval.m_arguments, retval.template nth_index<0>().end());
				} else {
					linear_int_helper<true>(retval);
				}
				return retval;
			}
		protected:
			/// Real power.
			/**
			* This method is written to work in conjunction with base_series::b_pow.
			*/
			template <class ArgsTuple>
			Derived real_pow(const double &y, const ArgsTuple &args_tuple) const {
				typedef typename Derived::term_type term_type;
				// Here we know that the cases of single coefficient, empty series and natural power have already been taken care of
				// in base_series::b_pow. We also know that if y is an integer, it must be -1.
				p_assert(!derived_const_cast->is_single_cf() and !derived_const_cast->empty());
				const int pow_n = (int)nearbyint(y);
				// If y is an integer, assert that it is -1 (i.e., that we are inverting).
				const bool is_inversion = std::abs(y - pow_n) < settings::numerical_zero();
				p_assert(!is_inversion or pow_n == -1);
				(void)is_inversion;
				term_type A(*derived_const_cast->template nth_index<0>().begin());
				// This is X, i.e., the original series without the leading term, which will then be divided by A.
				Derived XoverA(*derived_const_cast);
				XoverA.template term_erase<0>(args_tuple, XoverA.template nth_index<0>().begin());
				// Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
				term_type tmp_term;
				tmp_term.m_cf = A.m_cf.pow(-1, args_tuple);
				tmp_term.m_key = A.m_key.pow(-1, args_tuple);
				Derived Ainv;
				Ainv.insert(tmp_term, args_tuple, Ainv.template nth_index<0>().end());
				// Now let's compute X/A.
				XoverA.mult_by(Ainv, args_tuple);
				// Get the expansion limit from the truncator.
				size_t n;
				try {
					n = Derived::multiplier_type::truncator_type::power_series_limit(XoverA, args_tuple);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Poisson series is unsuitable for real exponentiation.\nThe reported error is: ")
									 + u.what());
				}
				return binomial_expansion(A, XoverA, y, n, args_tuple);
			}
		private:
			template <bool Flavour>
			void linear_int_helper(Derived &retval) const {
				// Handle the case in which the Poisson series is logically equivalent to a single polynomial with
				// linearly dependent integer coefficients.
				typedef typename Derived::term_type term_type;
				typedef typename term_type::cf_type cf_type;
				typedef typename term_type::key_type key_type;
				try {
					if (!derived_const_cast->is_single_cf()) {
						throw unsuitable("Series is not a linear combination of arguments with integer coefficients.");
					}
					// The size of the integer vector shall be the same as the poly arguments set's.
					std::vector<int> v(derived_const_cast->m_arguments.template get<0>().size());
					derived_const_cast->template nth_index<0>().begin()->m_cf.get_int_linear_combination(v);
					// Assign poly args of this as trig args of return value.
					retval.m_arguments.template get<1>() = derived_const_cast->m_arguments.template get<0>();
					term_type term;
					term.m_cf = cf_type((max_fast_int)1, retval.m_arguments);
					term.m_key.assign_int_vector(v);
					term.m_key.flavour() = Flavour;
					retval.insert(term, retval.m_arguments, retval.template nth_index<0>().end());
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("This Poisson series is not suitable for the application of circular functions. ") +
									 u.what());
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
