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

#ifndef PIRANHA_COMMON_POLYNOMIAL_CF_TOOLBOX_H
#define PIRANHA_COMMON_POLYNOMIAL_CF_TOOLBOX_H

#include <utility>
#include <vector>

#include "../integer_typedefs.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class common_polynomial_cf_toolbox
	{
		public:
			/// Return a single coefficient and a vector of integers representing the polynomial.
			template <class Cf>
			void get_int_linear_combination(std::pair<std::vector<Cf>, std::vector<max_fast_int> > &res) const {
				typedef typename Derived::const_sorted_iterator const_sorted_iterator;
				const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
				for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it) {
					const max_fast_int pos = it->m_key.linear_arg_position();
					if (pos >= 0) {
						p_assert(pos < (max_fast_int)res.second.size());
						res.second[(size_t)pos] = it->m_cf.get_int();
					} else {
						res.first.push_back(it->m_cf);
					}
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
