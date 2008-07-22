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

#ifndef PIRANHA_BASE_SERIES_PROBE_H
#define PIRANHA_BASE_SERIES_PROBE_H

namespace piranha
{
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool base_series<__PIRANHA_BASE_SERIES_TP>::is_single_cf() const
	{
		return (length() == 1 && begin()->m_key.is_unity());
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline double base_series<__PIRANHA_BASE_SERIES_TP>::norm(const ArgsTuple &args_tuple) const
	{
		const typename const_iterator::type it_f = end();
		double retval = 0;
		for (typename const_iterator::type it = begin(); it != it_f; ++it) {
			retval += it->m_cf.norm(args_tuple) * it->m_key.norm(args_tuple);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline typename base_series<__PIRANHA_BASE_SERIES_TP>::eval_type
	base_series<__PIRANHA_BASE_SERIES_TP>::eval(const double &t, const ArgsTuple &args_tuple) const
	{
		const typename const_iterator::type it_f = end();
		eval_type retval(0);
		for (typename const_iterator::type it = begin(); it != it_f; ++it) {
			eval_type tmp(it->m_cf.eval(t, args_tuple));
			tmp *= it->m_key.eval(t, args_tuple);
			retval += tmp;
		}
		return retval;
	}

	/// Return the number of elements of the series.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline size_t base_series<__PIRANHA_BASE_SERIES_TP>::length() const
	{
		return m_container.size();
	}

	/// Is series empty?
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline bool base_series<__PIRANHA_BASE_SERIES_TP>::empty() const
	{
		return m_container.empty();
	}

	/// Number of atoms in the series.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline size_t base_series<__PIRANHA_BASE_SERIES_TP>::atoms() const
	{
		size_t retval = 0;
		const typename const_iterator::type it_f = end();
		for (typename const_iterator::type it = begin(); it != it_f; ++it) {
			retval += it->m_cf.atoms() + it->m_key.atoms();
		}
		return retval;
	}
}

#endif
