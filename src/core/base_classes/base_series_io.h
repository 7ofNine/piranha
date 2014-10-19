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
#include <boost/numeric/conversion/cast.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../Psym.h"
#include "base_series_def.h"
#include "base_series_mp.h"
#include "base_series_tag.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Construct series from numerical quantity.
	/**
	 * Equivalent to creating a series from a coefficient with baseSeriesFromCf() and coefficient type term_type::cf_type, the coefficient being
	 * constructed with the provided number.
	 *
	 * @param[in] x number used for series construction.
	 * @param[in] argsTuple arguments tuple of the series.
	 *
	 * @return derived series constructed from given number.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSeriesFromNumber(const Number &x, const ArgsTuple &argsTuple)
	{
		return baseSeriesFromCf(typename TermType::cf_type(x, argsTuple), argsTuple);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::printTermsPlain(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		const const_iterator itf = end();
		const_iterator it = begin();
		while (it != itf) 
        {
			it->print_plain(stream, argsTuple);
			++it;
			if (it != itf) 
            {
				stream << separator;
			}
		}
	}


	// TODO: rework and fix the printing functions.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::genericPrintTermsPretty(std::ostream &stream, const Iterator &start, const Iterator &end,
		const ArgsTuple &argsTuple) const
	{
		const std::size_t max_length = settings::get_max_pretty_print_size();
		std::size_t count = 0;
		for (Iterator it = start; it != end; ++it) 
        {
			std::ostringstream tmp_stream;
			FromIterator<Iterator>::get(it)->print_pretty(tmp_stream,argsTuple);
			std::string tmp(tmp_stream.str());
			// If this is not the first term, we need to add the "+" sign if appropriate.
			if (it != start && !tmp.empty() && tmp[0] != '-') 
            {
				tmp.insert(tmp.begin(),'+');
			}
			count += tmp.size();
			if (count > max_length) 
            {
				std::for_each(tmp.begin(),
					tmp.begin() + boost::numeric_cast<std::string::iterator::difference_type>(max_length - (count - boost::numeric_cast<std::size_t>(tmp.size()))),
					stream << boost::lambda::_1);
				stream << "...";
				break;
			}

			if (!tmp.empty()) 
            {
				stream << tmp;
			}
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::printTermsPretty(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		if (empty()) 
        {
			stream << '0';
		} else 
        {
			try {
				const std::vector<typename Derived::TermType const *> s(derived_const_cast->template get_sorted_series<Derived>(argsTuple));
				genericPrintTermsPretty(stream, s.begin(), s.end(), argsTuple);

			} catch (const value_error &) {
				genericPrintTermsPretty(stream, begin(), end(), argsTuple);
			}
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::genericPrintTermsTEX(std::ostream &stream, const Iterator &start, const Iterator &end,
		const ArgsTuple &argsTuple) const
	{

		for (Iterator it = start; it != end; ++it) 
        {
			std::ostringstream tmp_stream;
			FromIterator<Iterator>::get(it)->print_tex(tmp_stream,argsTuple);
			std::string tmp(tmp_stream.str());
			// If this is not the first term, we need to add the "+" sign if appropriate.
			if (it != start && !tmp.empty() && tmp[0] != '-') 
            {
				tmp.insert(tmp.begin(),'+');
			}

			if (!tmp.empty()) 
            {
				stream << tmp;
			}
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::printTermsTEX(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		if (empty()) 
        {
			stream << '0';
		} else 
        {
			try {
				const std::vector<typename Derived::TermType const *> s(derived_const_cast->template get_sorted_series<Derived>(argsTuple));
				genericPrintTermsTEX(stream, s.begin(), s.end(), argsTuple);

			} catch (const value_error &)
            {
				genericPrintTermsTEX(stream, begin(), end(), argsTuple);
			}
		}
	}


	/// Constructor from Psym and from position in the arguments set.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseConstructFromPsym(const Psym &p, const int n,
			const ArgsTuple &argsTuple)
	{
		PIRANHA_ASSERT(derived_cast->empty());
		insert(TermType(typename TermType::cf_type(p, n, argsTuple), typename TermType::key_type(p, n, argsTuple)), argsTuple);
	}


	/// Begin of the series.
	/**
	 * @return iterator to the first term of the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<__PIRANHA_BASE_SERIES_TP>::begin() const
	{
		return m_container.begin();
	}

	/// End of the series.
	/**
	 * @return iterator to the last term of the series.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<__PIRANHA_BASE_SERIES_TP>::end() const
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
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSeriesFromKey(const Key &key, const ArgsTuple &argsTuple)
	{
		Derived retval;
		SeriesFromKeyImpl<Key, typename TermType::key_type>::run(retval, key, argsTuple);
		return retval;
	}


	/// Construct series from a cf.
	/**
	 * The returning series will consist of a single term with provided coefficient.
	 * If the provided coefficient type does not appear in the series, a compile-time error will occur.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Cf, class ArgsTuple>
	inline Derived BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSeriesFromCf(const Cf &cf, const ArgsTuple &argsTuple)
	{
		Derived retval;
		SeriesFromCfImpl<Cf, typename TermType::cf_type>::run(retval, cf, argsTuple);
		return retval;
	}


	/// Trivial destructor.
	/**
	 * No side effects. Contains compile-time checks.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline BaseSeries<__PIRANHA_BASE_SERIES_TP>::~BaseSeries()
	{
		PIRANHA_STATIC_CHECK((boost::is_base_of<BaseSeriesTag,Derived>::value), "Final series class must derive from BaseSeries class.");
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
