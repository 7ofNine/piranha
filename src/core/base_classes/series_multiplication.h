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

#define derived_const_cast static_cast<Derived const *>(this) // using the this of the derived class. It has to 
#define derived_cast static_cast<Derived *>(this)

#include <boost/type_traits/is_same.hpp>
#include <cstddef>
#include <iostream>

#include "../config.h"
#include "../settings.h"

// 
// for a series/polynomial to provide multiplication with another series even of its own type the 
// implementing class has to inherit from series_multiplication. BaseSeries::baseMultiply depends on it.
// i.e. why is it in a seperate class? Because it requires a Truncator? which is not necessary for a BaseSeries but BaseSeries wouldn't work without it. 
// At least the NullTruncator is required. Set it up as default truncator??
//
namespace piranha
{
	template <class Derived, class Multiplier, class Truncator>
	class series_multiplication
	{
		public:

			template <class ArgsTuple>
			int psIterations(int const start, int const step, ArgsTuple const &argsTuple) const
			{
                // the number of iterations needed depends on the truncation level
				return Multiplier::template get_type<Derived, Derived, ArgsTuple, Truncator>::truncator_type::powerSeriesIterations(*derived_const_cast, start, step, argsTuple);
			}


			template <class Series, class ArgsTuple>
			std::vector<typename Series::TermType const *> get_sorted_series(const ArgsTuple &argsTuple) const
			{
				static const bool check = boost::is_same<Series, Derived>::value;
                static_assert(check, "");

				return Multiplier::template get_type<Derived, Derived, ArgsTuple, Truncator>
                                            ::truncator_type::template getSortedPointerVector<Series, ArgsTuple>(*derived_const_cast, argsTuple);
			}


			// Multiply term-by-term with another series, and place the result into retval.
			// Preconditions:
			// - argsTuple must be the result of a merging of arguments between the two series being multiplied,
            // This is done e.g. in the named_series operator*=()
            //
			template <class Derived2, class ArgsTuple>
			void multiply_by_series(const Derived2 &s2, const ArgsTuple &argsTuple)
			{
				typedef typename Derived::const_iterator const_iterator;
				typedef typename Derived::TermType term_type;

				PIRANHA_DEBUG(std::cout << "Input lengths for series multiplication: " << derived_const_cast->length() << ','
					<< s2.length() << std::endl);

				// Don't do anything if this is empty. The "0" remains "0".
				if (derived_const_cast->empty()) 
				{
					return;
				}

				// If the other series is empty, clear the container and return. This represents a "0"
				if (s2.empty())
				{
					derived_cast->clearTerms();
					return;
				}


				const settings::MultiplicationAlgorithm algo = settings::getMultiplicationAlgorithm();
				// Optimize the cases of single coefficient series.
				if (s2.isSingleCf() && algo == settings::ALGORITHM_AUTOMATIC)
				{
					derived_cast->baseMultBy(s2.begin()->cf, argsTuple);
				} else if (derived_const_cast->isSingleCf() && algo == settings::ALGORITHM_AUTOMATIC)
				{
					Derived tmp;
					tmp.insertRange(s2.begin(), s2.end(), argsTuple);
					tmp.baseMultBy(derived_const_cast->begin()->cf, argsTuple);
					derived_cast->baseSwap(tmp);
				} else
				{   // the normal case. two series with several terms
					Derived retval;
					typename Multiplier::template get_type<Derived, Derived2, ArgsTuple, Truncator> multiplier(*derived_const_cast, s2, retval, argsTuple);
					multiplier.performMultiplication();
					derived_cast->baseSwap(retval);
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
