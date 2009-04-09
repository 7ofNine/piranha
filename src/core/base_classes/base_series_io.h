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

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "../settings.h"

namespace piranha
{
	/// Construct from numerical quantity.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::construct_from_number(const Number &x, const ArgsTuple &args_tuple)
	{
		// Make sure we are being called from an empty series.
		p_assert(empty());
		term_type term;
		term.m_cf = cf_type(x, args_tuple);
		insert(term, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::print_terms_plain(std::ostream &stream,
			const ArgsTuple &args_tuple, int limit) const
	{
		settings::setup_stream(stream);
		size_t j = 0, lim;
		if (limit < 0) {
			lim = length();
		} else {
			lim = static_cast<size_t>(limit);
		}
		const const_iterator it_f = end();
		for (const_iterator it = begin();it != it_f;++it) {
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
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::print_terms_pretty(std::ostream &stream,
			const ArgsTuple &args_tuple, int) const
	{
		using namespace boost::lambda;
		settings::setup_stream(stream);
		if (empty()) {
			stream << '0';
		} else {
			const size_t max_length = settings::get_max_pretty_print_size();
			size_t count = 0;
			const const_iterator it_f = end(), it_i = begin();
			for (const_iterator it = it_i; it != it_f; ++it) {
				std::ostringstream tmp_stream;
				settings::setup_stream(tmp_stream);
				it->print_pretty(tmp_stream,args_tuple);
				std::string tmp(tmp_stream.str());
				// If this is not the first term, we need to add the "+" sign if appropriate.
				if (it != it_i && !tmp.empty() && tmp[0] != '-') {
					tmp.insert(tmp.begin(),'+');
				}
				count += tmp.size();
				if (count > max_length) {
					std::for_each(tmp.begin(), tmp.begin() + max_length - (count - tmp.size()), stream << (boost::lambda::_1));
					stream << "...";
					break;
				}
				if (!tmp.empty()) {
					stream << tmp;
				}
			}
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::print_terms_latex(std::ostream &stream,
			const ArgsTuple &args_tuple, int limit) const
	{
// TODO: to be implemented.
	}

	/// Constructor from psym and from position in the arguments set.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::construct_from_psym_p(const psym_p &p, const int &n,
			const ArgsTuple &args_tuple)
	{
		p_assert(derived_cast->empty());
		insert(term_type(cf_type(p, n, args_tuple), key_type(p, n, args_tuple)), args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::iterator
	toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::begin()
	{
		return m_container.begin();
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::const_iterator
	toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::begin() const
	{
		return m_container.begin();
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::iterator
	toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::end()
	{
		return m_container.end();
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::const_iterator
	toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::end() const
	{
		return m_container.end();
	}
}

#endif
