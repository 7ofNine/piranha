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

#define DerivedConstCast static_cast<Derived const *>(this)
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
        static_assert(N >= 0, "Invalid arguments position in base Fourier series toolbox.");

		//protected:
		public:

			// Integrate supposing that the symbol is present in the fourier series..
			template <typename PosTuple, typename ArgsTuple>
			Derived baseIntegrate(const PosTuple &posTuple, const ArgsTuple &argsTuple) const
			{
                static_assert(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
					                 "Size mismatch between args tuple and pos tuple in Fourier series integration.");

				typedef typename Derived::const_iterator const_iterator;
				// Make sure that the position tuple contains just one symbol in position N and that
				// the symbol is actually present.
				PIRANHA_ASSERT(posTuple.template get<N>().size() == 1);
				PIRANHA_ASSERT(posTuple.template get<N>()[0].first);
				Derived retval;
				const std::size_t    pos  = posTuple.template get<N>()[0].second;
				const const_iterator  itf  = DerivedConstCast->end();

				for (const_iterator it = DerivedConstCast->begin(); it != itf; ++it) 
                {
					if (it->key[pos] == 0) 
                    {
						PIRANHA_THROW(value_error,"cannot procede with integration, one of the terms of the"
							          " Fourier series does not contain the symbol in its trigonometric arguments.");
					}

					typename Derived::TermType tmp(*it);
					if (it->key.getFlavour()) 
                    {
						tmp.cf.divideBy(  it->key[pos], argsTuple);

					} else 
                    {
						tmp.cf.divideBy( -it->key[pos], argsTuple);
					}

					tmp.key.setFlavour(!tmp.key.getFlavour());
					retval.insert(tmp, argsTuple);
				}

				return retval;
			}
	};
}

#undef DerivedConstCast
#undef derived_cast

#endif
