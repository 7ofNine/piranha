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
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../psym.h"
#include "../settings.h"
#include "base_series_def.h"
#include "series_factory.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Construct from numerical quantity.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::construct_from_number(const Number &x, const ArgsTuple &args_tuple)
	{
		// Make sure we are being called from an empty series.
		piranha_assert(empty());
		Derived tmp(base_series_from_cf(typename term_type::cf_type(x,args_tuple),args_tuple));
		base_swap(tmp);
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
					std::for_each(tmp.begin(), tmp.begin() + max_length - (count - tmp.size()), stream << boost::lambda::_1);
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
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::print_terms_tex(std::ostream &stream,
			const ArgsTuple &args_tuple, int) const
	{
		settings::setup_stream(stream);
		if (empty()) {
			stream << '0';
		} else {
			const const_iterator it_f = end(), it_i = begin();
			for (const_iterator it = it_i; it != it_f; ++it) {
				std::ostringstream tmp_stream;
				settings::setup_stream(tmp_stream);
				it->print_tex(tmp_stream,args_tuple);
				std::string tmp(tmp_stream.str());
				// If this is not the first term, we need to add the "+" sign if appropriate.
				if (it != it_i && !tmp.empty() && tmp[0] != '-') {
					tmp.insert(tmp.begin(),'+');
				}
				if (!tmp.empty()) {
					stream << tmp;
				}
			}
		}
	}

	/// Constructor from psym and from position in the arguments set.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_construct_from_psym(const psym &p, const int &n,
			const ArgsTuple &args_tuple)
	{
		piranha_assert(derived_cast->empty());
		insert(term_type(typename term_type::cf_type(p, n, args_tuple), typename term_type::key_type(p, n, args_tuple)), args_tuple);
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

	/// Construct series from a key.
	/**
	 * The returning series will consist of a single term with unitary numerical coefficient and key in its position.
	 * For instance, if the key is \f$ \cos\left( x \right) \f$ one in a Fourier series, the resulting series will be  \f$ 1 \cdot \cos\left( x \right) \f$;
	 * the exponent key \f$ x^2 y \f$ used to build a Poisson series, instead, will result instead in the series \f$ 1 \cdot x^2 y \cdot \cos\left( 0 \right) \f$.
	 * If the provided key type does not appear in the series, a compile-time error will occur.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Key, class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_series_from_key(const Key &key, const ArgsTuple &args_tuple)
	{
		Derived retval;
		series_from_key_impl<Key, typename term_type::key_type>::run(retval,key,args_tuple);
		return retval;
	}

	/// Construct series from a cf.
	/**
	 * The returning series will consist of a single term with provided coefficient
	 * If the provided coefficient type does not appear in the series, a compile-time error will occur.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Cf, class ArgsTuple>
	inline Derived toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_series_from_cf(const Cf &cf, const ArgsTuple &args_tuple)
	{
		Derived retval;
		series_from_cf_impl<Cf, typename term_type::cf_type>::run(retval,cf,args_tuple);
		return retval;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
