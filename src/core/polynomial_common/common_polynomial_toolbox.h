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

#ifndef PIRANHA_COMMON_POLYNOMIAL_TOOLBOX_H
#define PIRANHA_COMMON_POLYNOMIAL_TOOLBOX_H

#include <cmath>

#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_comparisons.h"
#include "../common_functors.h"
#include "../int_power_cache.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class common_polynomial_toolbox:
		public binomial_exponentiation_toolbox<Derived,term_key_degree_comparison>
	{
		public:
			// NOTICE: this will become a plain cache (i.e., not an EBS one) once
			// multinomial exponentiation for polynomials is implemented.
			template <class ArgsTuple>
			class sub_cache: public int_power_cache<Derived, base_series_arithmetics<Derived,ArgsTuple> >
			{
					typedef int_power_cache<Derived,
						base_series_arithmetics<Derived,ArgsTuple> > ancestor;
				public:
					sub_cache():ancestor::int_power_cache() {}
					void setup(const Derived &s, const ArgsTuple *args_tuple) {
						this->m_arith_functor.m_args_tuple = args_tuple;
						this->m_container[0] = Derived(static_cast<max_fast_int>(1),*args_tuple);
						this->m_container[1] = s;
					}
			};
			template <class ArgsTuple>
			double norm(const ArgsTuple &args_tuple) const {
				return std::abs(derived_const_cast->eval(0,args_tuple));
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
