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
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::pad_right(const ArgsTuple &args_tuple)
	{
		typedef typename Derived::term_type term_type;
		typedef typename Derived::const_iterator const_iterator;
		if (derived_const_cast->empty()) {
			return;
		}
		Derived retval;
		const const_iterator it_f = derived_const_cast->end();
		for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
        {
			term_type term(*it);
			term.m_cf.pad_right(args_tuple);
			term.m_key.pad_right(args_tuple);
			retval.insert(term, args_tuple);
		}

		derived_cast->base_swap(retval);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::swap(Derived &s2)
	{
		derived_cast->base_swap(s2);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class Layout, class ArgsTuple>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::apply_layout(const Layout &l, const ArgsTuple &args_tuple)
	{
		Derived retval;
		derived_cast->apply_layout_to_terms(l, retval, args_tuple);
		derived_cast->base_swap(retval);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::trim_test(TrimFlags &tf) const
	{
		derived_const_cast->trim_test_terms(tf);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline Derived cf_series<__PIRANHA_CF_SERIES_TP>::trim(const TrimFlags &tf, const ArgsTuple &args_tuple) const
	{
		Derived retval;
		derived_const_cast->trim_terms(tf, retval, args_tuple);
		return retval;
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
	inline RetSeries cf_series<__PIRANHA_CF_SERIES_TP>::sub(const PosTuple &p,
		SubCaches &s, const ArgsTuple &a) const
	{
		return derived_const_cast->template base_sub<RetSeries,typename Derived::sub_functor>(p, s, a);
	}


	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class Series, class ArgsTuple>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::split(std::vector<std::vector<Series> > &retval, const int &n, const ArgsTuple &args_tuple) const
	{
		derived_const_cast->base_split(retval, n, args_tuple);
	}
}

#endif
