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
#include <cstddef>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "../config.h"
#include "../exceptions.h"

#define DerivedConstCast static_cast<Derived const *>(this)

namespace piranha
{
	/// Base polynomial toolbox.
	/**
	 * N represents the position of polynomial arguments in the arguments tuple.
	 */
	template <int N, class Derived>
	class BasePolynomial
	{
			PIRANHA_STATIC_CHECK(N >= 0, "Invalid arguments position in base polynomial toolbox.");
		//protected:
		public:

			// Integrate supposing that the symbol is present in the polynomial.
			template <typename PositionTuple, typename ArgsTuple>
			Derived baseIntegrate(PositionTuple const &positionTuple, ArgsTuple const &argsTuple) const
			{
				PIRANHA_STATIC_CHECK(boost::tuples::length<PositionTuple>::value == boost::tuples::length<ArgsTuple>::value,
					"Size mismatch between args tuple and pos tuple in polynomial integration.");
			
				typedef typename Derived::const_iterator Iterator;
				typedef typename Derived::TermType::KeyType::DegreeType DegreeType;
				// Make sure that the position tuple contains just one symbol in position N and that
				// the symbol is actually present.
				PIRANHA_ASSERT(positionTuple.template get<N>().size() == 1);  // only one symbol we integrate over
				PIRANHA_ASSERT(positionTuple.template get<N>()[0].first);     // and it is present
				
                Derived retval;
				
                const std::size_t position = positionTuple.template get<N>()[0].second; // the position of the symbol in the argument/key 
				const Iterator       itEnd = DerivedConstCast->end();

				std::vector<DegreeType> tmpExponents(argsTuple.template get<N>().size());

				for (Iterator it = DerivedConstCast->begin(); it != itEnd; ++it) 
				{
					if (it->key[position] == -1) 
					{
						PIRANHA_THROW(value_error, "Exponent is -1 in integrand polynomial, cannot proceed"); // this would lead to logarithm
					}

					std::copy(it->key.begin(), it->key.end(), tmpExponents.begin());
					tmpExponents[position] += 1; // increment the exponent 

					typename Derived::TermType tmp(*it);

					tmp.key.resize(boost::numeric_cast<typename Derived::TermType::KeyType::size_type>(tmpExponents.size()));
					std::copy(tmpExponents.begin(), tmpExponents.end(), tmp.key.begin());
					
                    tmp.cf.divideBy(it->key[position] + 1, argsTuple);   // divide by n + 1
					
                    retval.insert(tmp, argsTuple);
				}

				return retval;
			}
	};
}

#undef DerivedConstCast

#endif
