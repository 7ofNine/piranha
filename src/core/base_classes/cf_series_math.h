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
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::invert_sign(const ArgsTuple &args_tuple)
	{
		// TODO: improve performance on this.
		typedef typename Derived::const_iterator const_iterator;
		typedef typename Derived::term_type term_type;
		Derived retval;
		const const_iterator it_f = derived_const_cast->end();
		for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) {
			term_type term(*it);
			term.m_cf.invert_sign(args_tuple);
			// No need to check, we are merging terms from this series.
			retval.template insert<false, true>(term, args_tuple);
		}
		derived_cast->base_swap(retval);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::add(const Derived &s, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_add(s,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::add(const double &x, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_add(x,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::add(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_add(q,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::add(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_add(z,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::subtract(const Derived &s, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_subtract(s,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::subtract(const double &x, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_subtract(x,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::subtract(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_subtract(q,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::subtract(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_subtract(z,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::mult_by(const Derived &s, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_mult_by(s,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::mult_by(const double &x, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_mult_by(x,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::mult_by(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_mult_by(q,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::mult_by(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_mult_by(z,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divide_by(const double &x, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_divide_by(x,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divide_by(const mp_rational &q, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_divide_by(q,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived &cf_series<__PIRANHA_CF_SERIES_TP>::divide_by(const mp_integer &z, const ArgsTuple &args_tuple)
	{
		return derived_cast->base_divide_by(z,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived cf_series<__PIRANHA_CF_SERIES_TP>::pow(const double &y, const ArgsTuple &args_tuple) const
	{
		return derived_const_cast->base_pow(y,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline Derived cf_series<__PIRANHA_CF_SERIES_TP>::pow(const mp_rational &q, const ArgsTuple &args_tuple) const
	{
		return derived_const_cast->base_pow(q,args_tuple);
	}

	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class Series, class PosTuple, class ArgsTuple>
	inline Series cf_series<__PIRANHA_CF_SERIES_TP>::partial(const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const
	{
		Series retval;
		Derived::base_partial(*derived_const_cast,retval,pos_tuple,args_tuple);
		return retval;
	}
}

#endif
