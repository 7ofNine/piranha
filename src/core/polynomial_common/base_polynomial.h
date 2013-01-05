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

#ifndef PIRANHA_BASE_POLYNOMIAL_H
#define PIRANHA_BASE_POLYNOMIAL_H

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <vector>

#include "../config.h"
#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Base polynomial toolbox.
	/**
	 * N represents the position of polynomial arguments in the arguments tuple.
	 */
	template <int N, class Derived>
	class BasePolynomial
	{
			p_static_check(N >= 0, "Invalid arguments position in base polynomial toolbox.");
		//protected:
		public:

			// Integrate supposing that the symbol is present in the polynomial.
			template <typename PosTuple, typename ArgsTuple>
			Derived base_integrate(const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const
			{
				p_static_check(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
					"Size mismatch between args tuple and pos tuple in polynomial integration.");
			
				typedef typename Derived::const_iterator const_iterator;
				typedef typename Derived::term_type::key_type::degree_type degree_type;
				// Make sure that the position tuple contains just one symbol in position N and that
				// the symbol is actually present.
				piranha_assert(pos_tuple.template get<N>().size() == 1);
				piranha_assert(pos_tuple.template get<N>()[0].first);
				Derived retval;
				const std::size_t pos = pos_tuple.template get<N>()[0].second;
				const const_iterator it_f = derived_const_cast->end();
				std::vector<degree_type> tmp_expos(args_tuple.template get<N>().size());

				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
				{
					if (it->m_key[pos] == -1) 
					{
						piranha_throw(value_error,"exponent is -1 in integrand polynomial, cannot proceed");
					}

					std::copy(it->m_key.begin(),it->m_key.end(),tmp_expos.begin());
					tmp_expos[pos] += 1;
					typename Derived::term_type tmp(*it);
					tmp.m_key.resize(boost::numeric_cast<typename Derived::term_type::key_type::size_type>(tmp_expos.size()));
					std::copy(tmp_expos.begin(),tmp_expos.end(),tmp.m_key.begin());
					tmp.m_cf.divide_by(it->m_key[pos] + 1,args_tuple);
					retval.insert(tmp,args_tuple);
				}

				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
