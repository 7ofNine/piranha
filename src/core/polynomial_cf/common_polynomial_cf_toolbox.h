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
	class common_polynomial_cf
	{
		public:
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct ei_sub_cache_selector {
				typedef typename Derived::TermType::CfType::
					template ei_sub_cache_selector<SubSeries, typename Derived::TermType::KeyType::
					template ei_sub_cache_selector<SubSeries, SubCachesCons, ArgsTuple>::type, ArgsTuple>::type type;
			};


			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &p, SubCaches &s, const ArgsTuple &argsTuple) const
			{
				return derived_const_cast->template baseSub<RetSeries, ei_sub_functor>(p, s, argsTuple);
			}


			template <class PosTuple, class ArgsTuple>
			Derived integrate(const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
			{
				return derived_const_cast->baseIntegrate(pos_tuple, argsTuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
