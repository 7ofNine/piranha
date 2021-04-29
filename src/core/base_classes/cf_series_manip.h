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

#ifndef PIRANHA_CF_SERIES_MANIP_H
#define PIRANHA_CF_SERIES_MANIP_H

#include <vector>

namespace piranha
{
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void CfSeries<PIRANHA_CF_SERIES_TP>::padRight(const ArgsTuple &argsTuple)
	{
		typedef typename Derived::TermType TermType;
		typedef typename Derived::const_iterator const_iterator;
		
		if (derived_const_cast->empty()) 
		{
			return;
		}

		Derived retval;
		const const_iterator itEnd = derived_const_cast->end();
		for (const_iterator it = derived_const_cast->begin(); it != itEnd; ++it) 
        {
			TermType term(*it);
			term.cf.padRight(argsTuple);
			term.key.padRight(argsTuple);

			retval.insert(term, argsTuple);
		}

		derived_cast->baseSwap(retval);
	}


	template <PIRANHA_CF_SERIES_TP_DECL>
	inline void CfSeries<PIRANHA_CF_SERIES_TP>::swap(Derived &series2)
	{
		derived_cast->baseSwap(series2);
	}

    //
    // apply a layout combining this with argsTuple
    //
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class LayoutTuple, class ArgsTuple>
	inline void CfSeries<PIRANHA_CF_SERIES_TP>::applyLayout(LayoutTuple const &layoutTuple, ArgsTuple const &argsTuple)
	{
		Derived retval;
		derived_cast->applyLayoutToTerms(layoutTuple, argsTuple, retval);
		derived_cast->baseSwap(retval);
	}

    //
    // test if symbol is used. If not it can be removed
    //
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void CfSeries<PIRANHA_CF_SERIES_TP>::trimTest(TrimFlags &trimFlags) const
	{
		derived_const_cast->trimTestTerms(trimFlags);
	}


    // 
    // remove unused symbols
    //
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline Derived CfSeries<PIRANHA_CF_SERIES_TP>::trim(const TrimFlags &trimFlags, const ArgsTuple &argsTuple) const
	{
		Derived retval;
		derived_const_cast->trimTerms(trimFlags, argsTuple, retval);
		return retval;
	}


	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
	inline RetSeries CfSeries<PIRANHA_CF_SERIES_TP>::sub(const PosTuple &positionTuple,
		SubCaches &subCaches, const ArgsTuple &argsTuple) const
	{
		return derived_const_cast->template baseSub<RetSeries, typename Derived::SubstitutionFunctor>(positionTuple, subCaches, argsTuple);
	}


	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class Series, class ArgsTuple>
	inline void CfSeries<PIRANHA_CF_SERIES_TP>::split(std::vector<std::vector<Series> > &retval, const int n, const ArgsTuple &argsTuple) const
	{
		derived_const_cast->baseSplit(retval, n, argsTuple);
	}
}

#endif
