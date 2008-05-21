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
					std::vector<max_fast_int> v(derived_const_cast->m_arguments.template get<0>().size());
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
