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

#include "../base_classes/cf_series.h"
#include "../common_functors.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"
#include "../polynomial_common/common_polynomial_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	// NOTE: this assumes that exponents are in position 0 of arguments tuple.
	class common_polynomial_cf_toolbox: public common_polynomial_toolbox<Derived>
	{
		public:
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct ei_sub_cache_selector {
				typedef typename Derived::term_type::cf_type::
					template ei_sub_cache_selector<SubSeries,typename Derived::term_type::key_type::
					template ei_sub_cache_selector<SubSeries,SubCachesCons,ArgsTuple>::type,ArgsTuple>::type type;
			};
			/// Return a single coefficient and a vector of integers representing the polynomial.
			template <int TargetPos, class Cf, class ArgsTuple>
			void get_int_linear_combination(std::pair<std::vector<Cf>, std::vector<max_fast_int> > &res,
											const ArgsTuple &args_tuple) const {
				typedef typename Derived::const_iterator const_iterator;
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) {
					const max_fast_int pos = it->m_key.linear_arg_position();
					if (pos >= 0) {
						// Let's find the corresponding symbol in the target vector of arguments.
						size_t i = 0;
						for (; i < args_tuple.template get<TargetPos>().size(); ++i) {
							if (args_tuple.template get<0>()[static_cast<size_t>(pos)] ==
									args_tuple.template get<TargetPos>()[i]) {
								break;
							}
						}
						p_assert(i != args_tuple.template get<TargetPos>().size());
						p_assert(i < res.second.size());
						res.second[i] = it->m_cf.get_int();
					} else {
						res.first.push_back(it->m_cf);
					}
				}
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &p, SubCaches &s, const ArgsTuple &args_tuple) const
			{
				return derived_const_cast->template base_sub<RetSeries,ei_sub_functor>(
					p,s,args_tuple
				);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
