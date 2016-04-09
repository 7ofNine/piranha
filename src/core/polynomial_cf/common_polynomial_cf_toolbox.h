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

#include <cstddef>
#include <utility>
#include <vector>

#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	// NOTE: this assumes that exponents are in position 0 of arguments tuple.
	template <class Derived>
	class CommonPolynomialCf
	{
		public:
			template <class SubstitutionSeries, class SubCachesCons, class ArgsTuple>
			class EiSubstitutionCacheSelector
            {
                public:

				typedef typename Derived::TermType::CfType::
					template EiSubstitutionCacheSelector<SubstitutionSeries, typename Derived::TermType::KeyType::
					template EiSubstitutionCacheSelector<SubstitutionSeries, SubCachesCons, ArgsTuple>::Type, ArgsTuple>::Type Type;
			};


			template <class RetSeries, class PositionTuple, class SubstitutionCaches, class ArgsTuple>
			RetSeries eiSubstitute(PositionTuple const &positionTuple, SubstitutionCaches &substitutionCaches, ArgsTuple const &argsTuple) const
			{
				return derived_const_cast->template baseSub<RetSeries, EiSubstitutionFunctor>(positionTuple, substitutionCaches, argsTuple);
			}


			template <class PositionTuple, class ArgsTuple>
			Derived integrate(PositionTuple const &positionTuple, ArgsTuple const &argsTuple) const
			{
				return derived_const_cast->baseIntegrate(positionTuple, argsTuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
