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

#include <boost/utility.hpp>
#include <cstddef>
#include <vector>

#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "../settings.h"
#include "named_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <bool Sign, class Derived2>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::mergeWithSeries(const Derived2 &s2)
	{
		// If we are merging with self, create a copy and call recursively.
		if ((void *)(boost::addressof(*derived_cast)) == (void *)(boost::addressof(s2)))
		{
			PIRANHA_DEBUG(std::cout << "Merging with self, performing a copy." << '\n');
			mergeWithSeries<Sign>(Derived(*derived_const_cast));

		} else 
		{
			mergeArgs(s2);
			if (Sign) 
			{
				derived_cast->baseAdd(s2, argumentsTuple);

			} else 
			{
				derived_cast->baseSubtract(s2, argumentsTuple);
			}
		}

		trim();
		return *derived_cast;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::multiplyBySeries(const Derived2 &s2)
	{
		// First we merge the arguments of the two series.
		mergeArgs(s2);
		// Then we perform the multiplication.
		derived_cast->baseMultBy(s2, argumentsTuple);
		trim();
		return *derived_cast;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator+=(const T &x)
	{
		return NamedSeriesAddSelector<T>::run(*derived_cast, x);
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator-=(const T &x)
	{
		return NamedSeriesSubtractSelector<T>::run(*derived_cast, x);
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator-() const
	{
		Derived retval(*derived_const_cast);
		retval *= -1;
		return retval;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::multiplyNumberHelper(const Number &x)
	{
		derived_cast->baseMultBy(x, argumentsTuple);
		trim();
		return *derived_cast;
	}

	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator*=(const T &x)
	{
		return NamedSeriesMultiplySelector<T>::run(*derived_cast, x);
	}

	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::divideNumberHelper(const Number &x)
	{
		derived_cast->baseDivideBy(x, argumentsTuple);
		trim();
		return *derived_cast;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline Derived & NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator/=(const T &x)
	{
		return divideNumberHelper(x);
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::pow(const double x) const
	{
		Derived retval(derived_const_cast->basePow(x, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim();
		return retval;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::pow(const mp_rational &q) const
	{
		Derived retval(derived_const_cast->basePow(q, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim();
		return retval;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::root(const int n) const
	{
		Derived retval(derived_const_cast->baseRoot(n, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim();
		return retval;
	}


	/// Partial derivative with respect to a piranha::Psym.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::partial(std::string const &name, int const n) const
	{
		typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type PositionTupleType;
		const Psym p(name);
		const PositionTupleType positionTuple = psyms2pos(VectorPsym(1, p), argumentsTuple);
		Derived retval(derived_const_cast->basePartial(n, positionTuple, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim();
		return retval;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
