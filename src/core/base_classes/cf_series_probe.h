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

#ifndef PIRANHA_CF_SERIES_PROBE_H
#define PIRANHA_CF_SERIES_PROBE_H

namespace piranha
{
    // check on insertability. What is the criterion? needs description
    // empty or if coefficient and key are insertable. determined by argsTuple?
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool CfSeries<PIRANHA_CF_SERIES_TP>::isInsertable(const ArgsTuple &argsTuple) const
	{
		return (derived_const_cast->empty() || (derived_const_cast->begin()->cf.isInsertable(argsTuple) &&
			                                    derived_const_cast->begin()->key.isInsertable(argsTuple)));
	}

    // check on padding needed. What is padding??
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool CfSeries<PIRANHA_CF_SERIES_TP>::needsPadding(const ArgsTuple &argsTuple) const
	{
		return (!derived_const_cast->empty() && (derived_const_cast->begin()->cf.needsPadding(argsTuple) ||
			                                     derived_const_cast->begin()->key.needsPadding(argsTuple)));
	}

    //empty series can be ignored
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool CfSeries<PIRANHA_CF_SERIES_TP>::isIgnorable(const ArgsTuple &) const
	{
		return (derived_const_cast->empty());
	}

    //norm of series
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline double CfSeries<PIRANHA_CF_SERIES_TP>::norm(const ArgsTuple &argsTuple) const
	{
		return (derived_const_cast->baseNorm(argsTuple));
	}

    //evaluate series for value t
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline typename TermEvalTypeDeterminer<Term>::Type CfSeries<PIRANHA_CF_SERIES_TP>::eval(const double t, const ArgsTuple &argsTuple) const
	{
		return (derived_const_cast->baseEval(t, argsTuple));
	}

    // equality check
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class T>
	inline bool CfSeries<PIRANHA_CF_SERIES_TP>::operator==(const T &x) const
	{
		return derived_const_cast->baseEqualTo(x);
	}

    //inequality check
	template <PIRANHA_CF_SERIES_TP_DECL>
	template <class T>
	inline bool CfSeries<PIRANHA_CF_SERIES_TP>::operator!=(const T &x) const
	{
		return !(operator==(x));
	}
}

#endif
