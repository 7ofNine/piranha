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
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool cf_series<__PIRANHA_CF_SERIES_TP>::is_insertable(const ArgsTuple &argsTuple) const
	{
		return (derived_const_cast->empty() || (derived_const_cast->begin()->cf.is_insertable(argsTuple) &&
			derived_const_cast->begin()->key.is_insertable(argsTuple)));
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool cf_series<__PIRANHA_CF_SERIES_TP>::needs_padding(const ArgsTuple &argsTuple) const
	{
		return (!derived_const_cast->empty() && (derived_const_cast->begin()->cf.needs_padding(argsTuple) ||
			derived_const_cast->begin()->key.needs_padding(argsTuple)));
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline bool cf_series<__PIRANHA_CF_SERIES_TP>::is_ignorable(const ArgsTuple &) const
	{
		return (derived_const_cast->empty());
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline double cf_series<__PIRANHA_CF_SERIES_TP>::norm(const ArgsTuple &argsTuple) const
	{
		return (derived_const_cast->base_norm(argsTuple));
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline typename TermEvalTypeDeterminer<Term>::type
		cf_series<__PIRANHA_CF_SERIES_TP>::eval(const double &t, const ArgsTuple &argsTuple) const
	{
		return (derived_const_cast->base_eval(t,argsTuple));
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class T>
	inline bool cf_series<__PIRANHA_CF_SERIES_TP>::operator==(const T &x) const
	{
		return derived_const_cast->base_equal_to(x);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class T>
	inline bool cf_series<__PIRANHA_CF_SERIES_TP>::operator!=(const T &x) const
	{
		return !(operator==(x));
	}
}

#endif
