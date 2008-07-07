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

#ifndef PIRANHA_BASE_SERIES_IO_H
#define PIRANHA_BASE_SERIES_IO_H

#include "../stream_manager.h"

namespace piranha
{
	/// Construct from numerical quantity.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::construct_from_number(const Number &x, const ArgsTuple &args_tuple)
	{
		// Make sure we are being called from an empty series.
		p_assert(derived_const_cast->template nth_index<0>().empty());
		term_type term;
		term.m_cf = cf_type(x, args_tuple);
		insert(term, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_plain(std::ostream &stream,
			const ArgsTuple &args_tuple, int limit) const
	{
		typedef typename Derived::const_sorted_iterator const_sorted_iterator;
		stream_manager::setup_print(stream);
		size_t j = 0, lim;
		if (limit < 0) {
			lim = derived_const_cast->template nth_index<0>().size();
		} else {
			lim = (size_t)limit;
		}
		const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
		for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin();it != it_f;++it) {
			if (j == lim) {
				break;
			}
			it->print_plain(stream, args_tuple);
			if (j < lim - 1) {
				stream << separator;
			}
			++j;
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_latex(std::ostream &stream,
			const ArgsTuple &args_tuple, int limit) const
	{
// TODO: to be implemented.
	}

	/// Constructor from psym and from position in the arguments set.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::construct_from_psym_p(const psym_p &p, const int &n,
			const ArgsTuple &args_tuple)
	{
		p_assert(derived_cast->template nth_index<0>().empty());
		insert(term_type(cf_type(p, n, args_tuple), key_type(p, n, args_tuple)), args_tuple);
	}
}

#endif
