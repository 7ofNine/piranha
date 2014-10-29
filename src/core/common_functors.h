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

#ifndef PIRANHA_COMMON_FUNCTORS_H
#define PIRANHA_COMMON_FUNCTORS_H

#include <algorithm>

#include "exceptions.h"

// NOTE: maybe better move these and the caches into own subdirectory.

namespace piranha
{
	// TODO: use std::pow instead of pow method and inv().
	template <class T>
	struct NamedSeriesArithmetics
	{
		T inv(const T &orig) const
		{
			return orig.pow(-1);
		}


		void multiply(T &orig, const T &other) const
		{
			orig *= other;
		}


		template <class U>
		T pow(const T &orig, const U &y) const
		{
			return orig.pow(y);
		}
	};
    

	template <class T, class ArgsTuple>
	struct BaseSeriesArithmetics {

		BaseSeriesArithmetics():argsTuple(0) {}
		
        T inv(const T &orig) const
		{
			PIRANHA_ASSERT(argsTuple);

			return orig.basePow(-1, *argsTuple);
		}

		void multiply(T &orig, const T &other) const
		{
			PIRANHA_ASSERT(argsTuple);

			orig.baseMultBy(other, *argsTuple);
		}

		template <class U>
		T pow(const T &orig, const U &y) const
		{
			PIRANHA_ASSERT(argsTuple);

			return orig.basePow(y, *argsTuple);
		}

		mutable ArgsTuple const *argsTuple;
	};


	struct EiSubFunctor {

		template <class RetSeries, class Element, class PosTuple, class SubCaches, class ArgsTuple>
		static RetSeries run(const Element &e, const PosTuple &posTuple, SubCaches &subCaches, const ArgsTuple &argsTuple)
        {
			return e.template eiSub<RetSeries>(posTuple, subCaches, argsTuple);
		}
	};
}

#endif
