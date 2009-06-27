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

#ifndef PIRANHA_SERIES_MULTIPLICATION_H
#define PIRANHA_SERIES_MULTIPLICATION_H

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_SERIES_MULTIPLICATION_TP_DECL class Derived, \
				class Multiplier, \
					class Truncator
#define __PIRANHA_SERIES_MULTIPLICATION_TP Derived,Multiplier,Truncator

#include <boost/type_traits/is_same.hpp>
#include <iostream>

#include "../config.h"
#include "../settings.h"
#include "toolbox.h"

namespace piranha
{
	template <__PIRANHA_SERIES_MULTIPLICATION_TP_DECL>
	struct series_multiplication {};

	template <__PIRANHA_SERIES_MULTIPLICATION_TP_DECL>
	class toolbox<series_multiplication<__PIRANHA_SERIES_MULTIPLICATION_TP> >
	{
		protected:
			template <class ArgsTuple>
			size_t psi_(const int &start, const int &step, const ArgsTuple &args_tuple) const
			{
				return Multiplier::template get_type<Derived, Derived, ArgsTuple,
					Truncator>::truncator_type::power_series_iterations(*derived_const_cast,
					start,step,args_tuple);
			}
			template <class Series, class ArgsTuple>
			std::vector<typename Series::term_type const *> get_sorted_series(const ArgsTuple &args_tuple) const
			{
				static const bool check = boost::is_same<Series,Derived>::value;
				p_static_check(check,"");
				return Multiplier::template get_type<Derived, Derived, ArgsTuple,
					Truncator>::truncator_type::template get_sorted_pointer_vector<Series,ArgsTuple>(*derived_const_cast,args_tuple);
			}
			// Multiply term-by-term with another series, and place the result into retval.
			// Preconditions:
			// - args_tuple must be the result of a merging of arguments between the two series being multiplied,
			template <class Derived2, class ArgsTuple>
			void multiply_by_series(const Derived2 &s2, const ArgsTuple &args_tuple)
			{
				typedef typename Derived::const_iterator const_iterator;
				typedef typename Derived::term_type term_type;
				__PDEBUG(std::cout << "Input lengths for series multiplication: " << derived_const_cast->length() << ','
					<< s2.length() << '\n');
				// Don't do anything if this is empty.
				if (derived_const_cast->empty()) {
					return;
				}
				// If the other series is empty, clear the container and return.
				if (s2.empty()) {
					derived_cast->clear_terms();
					return;
				}
				// Optimize the cases of single coefficient series.
				if (s2.is_single_cf()) {
					derived_cast->multiply_coefficients_by(s2.begin()->m_cf, args_tuple);
				} else if (derived_const_cast->is_single_cf()) {
					Derived tmp;
					tmp.insert_range(s2.begin(),s2.end(),args_tuple);
					tmp.multiply_coefficients_by(derived_const_cast->begin()->m_cf, args_tuple);
					derived_cast->base_swap(tmp);
				} else {
					Derived retval;
					typename Multiplier::template get_type<Derived, Derived2, ArgsTuple, Truncator>
						m(*derived_const_cast, s2, retval, args_tuple);
					m.perform_multiplication();
					derived_cast->base_swap(retval);
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast
#undef __PIRANHA_SERIES_MULTIPLICATION_TP_DECL
#undef __PIRANHA_SERIES_MULTIPLICATION_TP

#endif
