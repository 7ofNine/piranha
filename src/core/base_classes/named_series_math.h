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

#ifndef PIRANHA_NAMED_SERIES_MATH_H
#define PIRANHA_NAMED_SERIES_MATH_H

//#include <boost/utility.hpp>

#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "../settings.h"
#include "named_series_def.h"

#include <cstddef>
#include <memory>
#include <vector>

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
    //
    // add or subtract two series (lower level)
    // Sign: true  - addition
    //       false - subtraction 
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <bool Sign, class Derived2>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::mergeWithSeries(Derived2 const &series2)
	{
		// If we are merging with self, create a copy and call recursively.
		if ((void *)(std::addressof(*derived_cast)) == (void *)(std::addressof(series2)))
		{
			PIRANHA_DEBUG(std::cout << "Merging with self, performing a copy." << '\n');
			mergeWithSeries<Sign>(Derived(*derived_const_cast)); // create a copy and merge the copy

		} else 
		{
			mergeArgs(series2);  // combine the argument tuples first

			if (Sign) 
			{
				derived_cast->baseAdd(series2, argumentsTuple);      // add

			} else 
			{
				derived_cast->baseSubtract(series2, argumentsTuple); //subtract
			}
		}

		trim(); // remove unused symbols

		return *derived_cast;
	}


    //
    // multiply with another series (lower level)
	//
    template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::multiplyBySeries(Derived2 const &series2)
	{
		mergeArgs(series2); // First we merge the arguments of the two series.
		
		derived_cast->baseMultBy(series2, argumentsTuple); // Then we perform the multiplication. (the terms)
		trim();  // remove unused symbols

		return *derived_cast;
	}


    //
    //  += 
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator+=(T const &x)
	{
		return NamedSeriesAddSelector<T>::run(*derived_cast, x); // select dependend on T
	}


    //
    //  -=
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator-=(T const &x)
	{
		return NamedSeriesSubtractSelector<T>::run(*derived_cast, x); // select dependend on T
	}


    //
    //  unitary -
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator-() const
	{
		Derived retval(*derived_const_cast);
		retval *= -1;  // multiply with 1
		return retval;
	}


    //
    // multiply with x and trim result
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::multiplyNumberHelper(Number const &x)
	{
		derived_cast->baseMultBy(x, argumentsTuple);
		trim();   // remove unused symbols

		return *derived_cast;
	}


    //
    //  *= 
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator*=(T const &x)
	{
		return NamedSeriesMultiplySelector<T>::run(*derived_cast, x); // select dependent on T
	}


    //
    // divide by x and trim result
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::divideNumberHelper(Number const &x)
	{
		derived_cast->baseDivideBy(x, argumentsTuple);
		trim();  // remove unused sybols

		return *derived_cast;
	}


    //
    //    /=
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator/=(T const &x)
	{
		return divideNumberHelper(x);
	}


    //
    //  power(double) of series and trim
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::pow(double const x) const
	{
		Derived retval(derived_const_cast->basePow(x, argumentsTuple)); //power
		retval.argumentsTuple = argumentsTuple;
		retval.trim();  // remove unused symbols

		return retval;
	}


    //
    //  power(rational) of series and trim
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::pow(mp_rational const &q) const
	{
		Derived retval(derived_const_cast->basePow(q, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim(); // remove unused symbols

		return retval;
	}


    //
    // n-th root of series and trim
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::root(int const n) const
	{
		Derived retval(derived_const_cast->baseRoot(n, argumentsTuple)); // root
		retval.argumentsTuple = argumentsTuple;
		retval.trim(); // remove unused symbols

		return retval;
	}


    //
	// n-th Partial derivative with respect to variable "name".
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::partial(std::string const &name, int const n) const
	{
		typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type PositionTupleType;
		// IS that correct. Doesn't Psym(Name) overwright the original settings if the name already exists??
		const PositionTupleType positionTuple = psyms2pos(VectorPsym(1, Psym(name)), argumentsTuple); // get position in argsTuple

		Derived retval(derived_const_cast->basePartial(n, positionTuple, argumentsTuple)); // partial derivative on terms
		retval.argumentsTuple = argumentsTuple;
		retval.trim(); // remove unused symbols

		return retval;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
