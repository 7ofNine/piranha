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

#ifndef PIRANHA_POWER_SERIES_H
#define PIRANHA_POWER_SERIES_H

#include <algorithm> // For max_element and min_element.
#include <boost/static_assert.hpp>

#include "../integer_typedefs.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Power series toolbox.
	template <int ExpoArgsPosition, int ExpoTermPosition, class Derived>
	class power_series
	{
			BOOST_STATIC_ASSERT(ExpoArgsPosition >= 0);
			template <class Term>
			class degree_binary_predicate
			{
				public:
					bool operator()(const Term &t1, const Term &t2) const {
						return (t1.template get<ExpoTermPosition>().degree() < t2.template get<ExpoTermPosition>().degree());
					}
			};
			template <class Term>
			class min_degree_binary_predicate
			{
				public:
					bool operator()(const Term &t1, const Term &t2) const {
						return (t1.template get<ExpoTermPosition>().min_degree() <
								t2.template get<ExpoTermPosition>().min_degree());
					}
			};
		public:
			static const int expo_args_position = ExpoArgsPosition;
			static const int expo_term_position = ExpoTermPosition;
			/// Get the degree of the power series.
			max_fast_int degree() const {
				if (derived_const_cast->empty()) {
					return 0;
				}
				const typename Derived::const_iterator result(std::max_element(
							derived_const_cast->begin(),
							derived_const_cast->end(),
							degree_binary_predicate<typename Derived::term_type>()
						));
				return result->template get<ExpoTermPosition>().degree();
			}
			/// Get the minimum degree of the power series.
			max_fast_int min_degree() const {
				if (derived_const_cast->empty()) {
					return 0;
				}
				const typename Derived::const_iterator result(std::min_element(
							derived_const_cast->begin(),
							derived_const_cast->end(),
							min_degree_binary_predicate<typename Derived::term_type>()
						));
				return result->template get<ExpoTermPosition>().min_degree();
			}
			void upload_min_exponents(std::vector<max_fast_int> &v) const {
				p_assert(!derived_const_cast->empty());
				typedef typename Derived::const_iterator const_iterator;
				const const_iterator it_f = derived_const_cast->end();
				const_iterator it = derived_const_cast->begin();
				it->template get<ExpoTermPosition>().upload_min_exponents(v);
				for (; it != it_f; ++it) {
					it->template get<ExpoTermPosition>().test_min_exponents(v);
				}
			}
			/// Get a vector containing the minimum exponents of the series.
			template <class ArgsTuple>
			std::vector<max_fast_int> min_exponents(const ArgsTuple &args_tuple) const {
				std::vector<max_fast_int> retval(args_tuple.template get<ExpoArgsPosition>().size());
				upload_min_exponents(retval);
				return retval;
			}
			template <class ArgsTuple>
			max_fast_int min_expo_of(const size_t &pos, const ArgsTuple &args_tuple) const {
				std::vector<max_fast_int> tmp(min_exponents(args_tuple));
				if (pos >= tmp.size()) {
					return 0;
				} else {
					return tmp[pos];
				}
			}
			void test_min_exponents(std::vector<max_fast_int> &v) const {
				typedef typename Derived::const_iterator const_iterator;
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) {
					it->template get<ExpoTermPosition>().test_min_exponents(v);
				}
			}
			// Return true if the minimum exponents are smaller than those specified in the limits vector.
			template <class ArgsTuple>
			bool test_expo_limits(const std::vector<std::pair<size_t, max_fast_int> > &v, const ArgsTuple &args_tuple) const {
				const std::vector<max_fast_int> min_expo(min_exponents(args_tuple));
				const size_t size = v.size();
				for (size_t i = 0; i < size; ++i) {
					p_assert(v[i].first < min_expo.size());
					if (v[i].second < min_expo[v[i].first]) {
						return false;
					}
				}
				return true;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
