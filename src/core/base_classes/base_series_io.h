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
#include <boost/type_traits/is_base_of.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../psym.h"
#include "base_series_def.h"
#include "base_series_mp.h"
#include "base_series_tag.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Construct series from numerical quantity.
	/**
	 * Equivalent to creating a series from a coefficient with base_series_from_cf() and coefficient type term_type::cf_type, the coefficient being
	 * constructed with the provided number.
	 *
	 * @param[in] x number used for series construction.
	 * @param[in] args_tuple arguments tuple of the series.
	 *
	 * @return derived series constructed from given number.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::base_series_from_number(const Number &x, const ArgsTuple &args_tuple)
	{
		return base_series_from_cf(typename term_type::cf_type(x,args_tuple),args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_plain(std::ostream &stream, const ArgsTuple &args_tuple) const
	{
		const const_iterator it_f = end();
		const_iterator it = begin();
		while (it != it_f) {
			it->print_plain(stream, args_tuple);
			++it;
			if (it != it_f) {
				stream << separator;
			}
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::generic_print_terms_pretty(std::ostream &stream, const Iterator &start, const Iterator &end,
		const ArgsTuple &args_tuple) const
	{
		const std::size_t max_length = settings::get_max_pretty_print_size();
		std::size_t count = 0;
		for (Iterator it = start; it != end; ++it) {
			std::ostringstream tmp_stream;
			it_getter<Iterator>::get(it)->print_pretty(tmp_stream,args_tuple);
			std::string tmp(tmp_stream.str());
			// If this is not the first term, we need to add the "+" sign if appropriate.
			if (it != start && !tmp.empty() && tmp[0] != '-') {
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

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_pretty(std::ostream &stream, const ArgsTuple &args_tuple) const
	{
		if (empty()) {
			stream << '0';
		} else {
			try {
				const std::vector<typename Derived::term_type const *> s(derived_const_cast->template get_sorted_series<Derived>(args_tuple));
				generic_print_terms_pretty(stream,&(*s.begin()),&(*s.end()),args_tuple);
			} catch (const value_error &) {
				generic_print_terms_pretty(stream,begin(),end(),args_tuple);
			}
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::generic_print_terms_tex(std::ostream &stream, const Iterator &start, const Iterator &end,
		const ArgsTuple &args_tuple) const
	{

		for (Iterator it = start; it != end; ++it) {
			std::ostringstream tmp_stream;
			it_getter<Iterator>::get(it)->print_tex(tmp_stream,args_tuple);
			std::string tmp(tmp_stream.str());
			// If this is not the first term, we need to add the "+" sign if appropriate.
			if (it != start && !tmp.empty() && tmp[0] != '-') {
				tmp.insert(tmp.begin(),'+');
			}
			if (!tmp.empty()) {
				stream << tmp;
			}
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_tex(std::ostream &stream, const ArgsTuple &args_tuple) const
	{
		if (empty()) {
			stream << '0';
		} else {
			try {
				const std::vector<typename Derived::term_type const *> s(derived_const_cast->template get_sorted_series<Derived>(args_tuple));
				generic_print_terms_tex(stream,&(*s.begin()),&(*s.end()),args_tuple);
			} catch (const value_error &) {
				generic_print_terms_tex(stream,begin(),end(),args_tuple);
			}
		}
	}

	/// Constructor from psym and from position in the arguments set.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::base_construct_from_psym(const psym &p, const int &n,
			const ArgsTuple &args_tuple)
	{
		piranha_assert(derived_cast->empty());
		insert(term_type(typename term_type::cf_type(p, n, args_tuple), typename term_type::key_type(p, n, args_tuple)), args_tuple);
	}

	/// Begin of the series.
	/**
	 * @return iterator to the first term of the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename base_series<__PIRANHA_BASE_SERIES_TP>::const_iterator
	base_series<__PIRANHA_BASE_SERIES_TP>::begin() const
	{
		return m_container.begin();
	}

	/// End of the series.
	/**
	 * @return iterator to the last term of the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename base_series<__PIRANHA_BASE_SERIES_TP>::const_iterator
	base_series<__PIRANHA_BASE_SERIES_TP>::end() const
	{
		return m_container.end();
	}

	/// Construct series from a key.
	/**
	 * The returning series will consist of a single term with unitary numerical coefficient and key in its position.
	 * For instance, if the key is \f$ \cos\left( x \right) \f$ in a Fourier series, the resulting series will be  \f$ 1 \cdot \cos\left( x \right) \f$;
	 * the exponent key \f$ x^2 y \f$ used to build a Poisson series, instead, will result instead in the series \f$ 1 \cdot x^2 y \cdot \cos\left( 0 \right) \f$.
	 * If the provided key type does not appear in the series, a compile-time error will occur.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Key, class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::base_series_from_key(const Key &key, const ArgsTuple &args_tuple)
	{
		Derived retval;
		series_from_key_impl<Key, typename term_type::key_type>::run(retval,key,args_tuple);
		return retval;
	}

	/// Construct series from a cf.
	/**
	 * The returning series will consist of a single term with provided coefficient.
	 * If the provided coefficient type does not appear in the series, a compile-time error will occur.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Cf, class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::base_series_from_cf(const Cf &cf, const ArgsTuple &args_tuple)
	{
		Derived retval;
		series_from_cf_impl<Cf, typename term_type::cf_type>::run(retval,cf,args_tuple);
		return retval;
	}

	/// Trivial destructor.
	/**
	 * No side effects. Contains compile-time checks.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline base_series<__PIRANHA_BASE_SERIES_TP>::~base_series()
	{
		p_static_check((boost::is_base_of<base_series_tag,Derived>::value),"Final series class must derive from base_series class.");
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
