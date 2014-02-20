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

#ifndef PIRANHA_BASE_FOURIER_SERIES_H
#define PIRANHA_BASE_FOURIER_SERIES_H

#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <vector>

#include "../config.h"
#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Base fourier series toolbox.
	/**
	 * N represents the position of trigonometric arguments in the arguments tuple.
	 */
	template <int N, class Derived>
	class BaseFourierSeries
	{
			p_static_check(N >= 0, "Invalid arguments position in base Fourier series toolbox.");

		//protected:
		public:

			// Integrate supposing that the symbol is present in the fourier series..
			template <typename PosTuple, typename ArgsTuple>
			Derived base_integrate(const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
			{
				p_static_check(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
					"Size mismatch between args tuple and pos tuple in Fourier series integration.");

				typedef typename Derived::const_iterator const_iterator;
				// Make sure that the position tuple contains just one symbol in position N and that
				// the symbol is actually present.
				piranha_assert(pos_tuple.template get<N>().size() == 1);
				piranha_assert(pos_tuple.template get<N>()[0].first);
				Derived retval;
				const std::size_t pos = pos_tuple.template get<N>()[0].second;
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
                {
					if (it->m_key[pos] == 0) 
                    {
						piranha_throw(value_error,"cannot procede with integration, one of the terms of the "
							" Fourier series does not contain the symbol in its trigonometric arguments.");
					}

					typename Derived::term_type tmp(*it);
					if (it->m_key.get_flavour()) 
                    {
						tmp.m_cf.divide_by(  it->m_key[pos], argsTuple);

					} else 
                    {
						tmp.m_cf.divide_by(- it->m_key[pos], argsTuple);
					}

					tmp.m_key.set_flavour(!tmp.m_key.get_flavour());
					retval.insert(tmp, argsTuple);
				}

				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
