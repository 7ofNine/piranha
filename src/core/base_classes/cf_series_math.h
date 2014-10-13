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

#ifndef PIRANHA_CF_SERIES_MATH_H
#define PIRANHA_CF_SERIES_MATH_H

#include "../math.h"

namespace piranha
{
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::invert_sign(const ArgsTuple &argsTuple)
	{
		// TODO: improve performance on this.
		typedef typename Derived::const_iterator const_iterator;
		typedef typename Derived::TermType term_type;
		Derived retval;
		const const_iterator it_f = derived_const_cast->end();
		for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it)
        {
			term_type term(*it);
			term.cf.invert_sign(argsTuple);
			// No need to check, we are merging terms from this series.
			retval.template insert<false, true>(term, argsTuple);
		}
		derived_cast->base_swap(retval);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::add(const T &x, const ArgsTuple &argsTuple)
	{
		return derived_cast->base_add(x, argsTuple);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::subtract(const T &x, const ArgsTuple &argsTuple)
	{
		return derived_cast->base_subtract(x, argsTuple);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::mult_by(const T &x, const ArgsTuple &argsTuple)
	{
		return derived_cast->base_mult_by(x, argsTuple);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class T, class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divideBy(const T &x, const ArgsTuple &argsTuple)
	{
		return derived_cast->base_divide_by(x, argsTuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived cf_series<__PIRANHA_CF_SERIES_TP>::pow(const double &y, const ArgsTuple &argsTuple) const
	{
		return derived_const_cast->base_pow(y, argsTuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived cf_series<__PIRANHA_CF_SERIES_TP>::pow(const mp_rational &q, const ArgsTuple &argsTuple) const
	{
		return derived_const_cast->base_pow(q, argsTuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class Series, class PosTuple, class ArgsTuple>
	inline Series cf_series<__PIRANHA_CF_SERIES_TP>::partial(const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
	{
		Series retval;
		Derived::base_partial(*derived_const_cast, retval, pos_tuple, argsTuple);
		return retval;
	}
}

#endif
