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

#ifndef PIRANHA_BASE_POWER_SERIES_H
#define PIRANHA_BASE_POWER_SERIES_H

#include <algorithm> // For max_element and min_element.

#include "../config.h"
#include "../exceptions.h"
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <int ExpoArgsPosition, int ExpoTermPosition, class Derived>
	struct base_power_series {};

	/// Power series toolbox.
	template <int ExpoArgsPosition, int ExpoTermPosition, class Derived>
	class toolbox<base_power_series<ExpoArgsPosition,ExpoTermPosition,Derived> >
	{
			p_static_check(ExpoArgsPosition >= 0, "Invalid expo args position.");
			template <class Term>
			struct degree_binary_predicate
			{
				bool operator()(const Term &t1, const Term &t2) const {
					return (t1.template get<ExpoTermPosition>().degree() < t2.template get<ExpoTermPosition>().degree());
				}
			};
			template <class Term>
			struct min_degree_binary_predicate
			{
				bool operator()(const Term &t1, const Term &t2) const {
					return (t1.template get<ExpoTermPosition>().min_degree() <
						t2.template get<ExpoTermPosition>().min_degree());
				}
			};
			template <class Term, class PosTuple>
			struct partial_degree_binary_predicate
			{
				partial_degree_binary_predicate(const PosTuple &p):m_p(p) {}
				bool operator()(const Term &t1, const Term &t2) const {
					return (t1.template get<ExpoTermPosition>().partial_degree(m_p) < t2.template get<ExpoTermPosition>().partial_degree(m_p));
				}
				const PosTuple &m_p;
			};
			template <class Term, class PosTuple>
			struct partial_min_degree_binary_predicate
			{
				partial_min_degree_binary_predicate(const PosTuple &p):m_p(p) {}
				bool operator()(const Term &t1, const Term &t2) const {
					return (t1.template get<ExpoTermPosition>().partial_min_degree(m_p) <
						t2.template get<ExpoTermPosition>().partial_min_degree(m_p));
				}
				const PosTuple &m_p;
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
		protected:
			/// Get the degree of the power series for specific variables.
			template <class PosTuple>
			int base_partial_degree(const PosTuple &pos_tuple) const {
				if (derived_const_cast->empty()) {
					return 0;
				}
				const typename Derived::const_iterator result(std::max_element(
							derived_const_cast->begin(),
							derived_const_cast->end(),
							partial_degree_binary_predicate<typename Derived::term_type,PosTuple>(pos_tuple)
						));
				return result->template get<ExpoTermPosition>().partial_degree(pos_tuple);
			}
			/// Get the mininum degree of the power series for specific variables.
			template <class PosTuple>
			int base_partial_min_degree(const PosTuple &pos_tuple) const {
				if (derived_const_cast->empty()) {
					return 0;
				}
				const typename Derived::const_iterator result(std::min_element(
							derived_const_cast->begin(),
							derived_const_cast->end(),
							partial_min_degree_binary_predicate<typename Derived::term_type,PosTuple>(pos_tuple)
						));
				return result->template get<ExpoTermPosition>().partial_min_degree(pos_tuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
