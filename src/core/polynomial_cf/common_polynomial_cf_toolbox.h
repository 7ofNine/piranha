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
#include "../polynomial_common/common_polynomial_toolbox.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	// NOTE: this assumes that exponents are in position 0 of arguments tuple.
	class common_polynomial_cf_toolbox: public common_polynomial_toolbox<Derived>
	{
			typedef typename cf_series<Derived>::template reference_proxy<Derived> proxy_ancestor;
		public:
			class proxy: public proxy_ancestor
			{
				public:
					proxy(const Derived &s): proxy_ancestor(s), m_min_expos_cached(false),m_norm_cached(false) {
						m_min_degree = proxy_ancestor::m_ptr->min_degree();
					}
					template <class ArgsTuple>
					piranha::max_fast_int min_expo_of(const size_t &n, const ArgsTuple &args_tuple) const {
						if (!m_min_expos_cached) {
							m_min_expos = proxy_ancestor::m_ptr->min_exponents(args_tuple);
							m_min_expos_cached = true;
						}
						if (n >= m_min_expos.size()) {
							return 0;
						} else {
							return m_min_expos[n];
						}
					}
					const piranha::max_fast_int &min_degree() const {
						return m_min_degree;
					}
					template <class ArgsTuple>
					const double &norm(const ArgsTuple &args_tuple) const {
						if (!m_norm_cached) {
							m_norm = proxy_ancestor::m_ptr->norm(args_tuple);
							m_norm_cached = true;
						}
						return m_norm;
					}
					template <int TargetPos, class Cf, class ArgsTuple>
					void get_int_linear_combination(std::pair<std::vector<Cf>, std::vector<max_fast_int> > &res,
											const ArgsTuple &args_tuple) const {
						proxy_ancestor::m_ptr->template get_int_linear_combination<TargetPos>(res,args_tuple);
					}
					template <class ArgsTuple>
					std::vector<max_fast_int> min_exponents(const ArgsTuple &args_tuple) const {
						return proxy_ancestor::m_ptr->min_exponents(args_tuple);
					}
					template <class ArgsTuple>
					Derived besselJ(const max_fast_int &order, const ArgsTuple &args_tuple) const {
						return proxy_ancestor::m_ptr->besselJ(order,args_tuple);
					}
					static const int expo_args_position = Derived::expo_args_position;
				private:
					mutable bool								m_min_expos_cached;
					mutable bool								m_norm_cached;
					mutable std::vector<piranha::max_fast_int>	m_min_expos;
					piranha::max_fast_int						m_min_degree;
					mutable double								m_norm;
			};
			/// Return a single coefficient and a vector of integers representing the polynomial.
			template <int TargetPos, class Cf, class ArgsTuple>
			void get_int_linear_combination(std::pair<std::vector<Cf>, std::vector<max_fast_int> > &res,
											const ArgsTuple &args_tuple) const {
				typedef typename Derived::const_iterator::type const_iterator;
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
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
