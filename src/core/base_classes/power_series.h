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

#include "../p_assert.h"
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <int ExpoArgsPosition, int ExpoTermPosition, class Derived>
	struct power_series {};

	/// Power series toolbox.
	template <int ExpoArgsPosition, int ExpoTermPosition, class Derived>
	class toolbox<power_series<ExpoArgsPosition,ExpoTermPosition,Derived> >
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
			int degree() const {
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
			int min_degree() const {
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
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
