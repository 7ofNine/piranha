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

#include <type_traits>

#define DCC static_cast<Derived const *>(this)
#define DC static_cast<Derived *>(this)

namespace piranha
{
	/// Construct series from numerical quantity.
	/**
	 * Equivalent to creating a series from a coefficient with baseSeriesFromCf() and coefficient type term_type::CfType, the coefficient being
	 * constructed with the provided number.
	 *
	 * @param[in] x number used for series construction.
	 * @param[in] argsTuple arguments tuple of the series.
	 *
	 * @return derived series constructed from given number.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Number, class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSeriesFromNumber(const Number &x, const ArgsTuple &argsTuple)
	{
		return baseSeriesFromCf(typename TermType::CfType(x, argsTuple), argsTuple);
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::printTermsPlain(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		const const_iterator itf = end();
		const_iterator       it  = begin();

		while (it != itf) 
        {
			it->printPlain(stream, argsTuple);
			++it;
			if (it != itf) 
            {
				stream << separator;
			}
		}
	}


	// TODO: rework and fix the printing functions.
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::genericPrintTermsPretty(std::ostream &stream, const Iterator &start, const Iterator &end, const ArgsTuple &argsTuple) const
	{
		const std::size_t maxLength = settings::get_max_pretty_print_size();
		std::size_t count = 0;
		for (Iterator it = start; it != end; ++it) 
        {
			std::ostringstream tmpStream;
			FromIterator<Iterator>::get(it)->printPretty(tmpStream, argsTuple);
			std::string tmp(tmpStream.str());
			// If this is not the first term, we need to add the "+" sign if appropriate.
			if (it != start && !tmp.empty() && tmp[0] != '-') 
            {
				tmp.insert(tmp.begin(), '+');
			}

			count += tmp.size();
			if (count > maxLength) 
            {
				std::for_each(tmp.begin(),
					tmp.begin() + boost::numeric_cast<std::string::iterator::difference_type>(maxLength - (count - boost::numeric_cast<std::size_t>(tmp.size()))),
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


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::printTermsPretty(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		if (empty()) 
        {
			stream << '0';

		} else 
        {
			try {
				const std::vector<typename Derived::TermType const *> s(DCC->template get_sorted_series<Derived>(argsTuple));
				genericPrintTermsPretty(stream, s.begin(), s.end(), argsTuple);

			} catch (const value_error &)
            {
				genericPrintTermsPretty(stream, begin(), end(), argsTuple);
			}
		}
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::genericPrintTermsTEX(std::ostream &stream, const Iterator &start, const Iterator &end, const ArgsTuple &argsTuple) const
	{
		for (Iterator it = start; it != end; ++it) 
        {
			std::ostringstream tmpStream;
			FromIterator<Iterator>::get(it)->printTex(tmpStream, argsTuple);
			std::string tmp(tmpStream.str());
			// If this is not the first term, we need to add the "+" sign if appropriate.
			if (it != start && !tmp.empty() && tmp[0] != '-') 
            {
				tmp.insert(tmp.begin(), '+');
			}

			if (!tmp.empty()) 
            {
				stream << tmp;
			}
		}
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::printTermsTEX(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		if (empty()) 
        {
			stream << '0';
		} else 
        {
			try {
				const std::vector<typename Derived::TermType const *> s(DCC->template get_sorted_series<Derived>(argsTuple));
				genericPrintTermsTEX(stream, s.begin(), s.end(), argsTuple);

			} catch (const value_error &)
            {
				genericPrintTermsTEX(stream, begin(), end(), argsTuple);
			}
		}
	}


	/// Constructor from Psym and from position in the arguments set.
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::baseConstructFromPsym(const Psym &p, const int n, const ArgsTuple &argsTuple)
	{
		PIRANHA_ASSERT(DC->empty());
		insert(TermType(typename TermType::CfType(p, n, argsTuple), typename TermType::KeyType(p, n, argsTuple)), argsTuple);
	}


	/// Begin of the series.
	/**
	 * @return iterator to the first term of the series.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<PIRANHA_BASE_SERIES_TP>::begin() const
	{
		return container.begin();
	}

	/// End of the series.
	/**
	 * @return iterator to the last term of the series.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<PIRANHA_BASE_SERIES_TP>::end() const
	{
		return container.end();
	}


	/// Construct series from a key.
	/**
	 * The returning series will consist of a single term with unitary numerical coefficient and key in its position.
	 * For instance, if the key is \f$ \cos\left( x \right) \f$ in a Fourier series, the resulting series will be  \f$ 1 \cdot \cos\left( x \right) \f$;
	 * the exponent key \f$ x^2 y \f$ used to build a Poisson series, instead, will result instead in the series \f$ 1 \cdot x^2 y \cdot \cos\left( 0 \right) \f$.
	 * If the provided key type does not appear in the series, a compile-time error will occur.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Key, class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSeriesFromKey(const Key &key, const ArgsTuple &argsTuple)
	{
		Derived retval;
		SeriesFromKeyImpl<Key, typename TermType::KeyType>::run(retval, key, argsTuple);
		return retval;
	}


	/// Construct series from a cf.
	/**
	 * The returning series will consist of a single term with provided coefficient.
	 * If the provided coefficient type does not appear in the series, a compile-time error will occur.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Cf, class ArgsTuple>
	inline Derived BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSeriesFromCf(const Cf &cf, const ArgsTuple &argsTuple)
	{
		Derived retval;
		SeriesFromCfImpl<Cf, typename TermType::CfType>::run(retval, cf, argsTuple);
		return retval;
	}


	// Trivial destructor.
	
	// No side effects. Contains compile-time checks.
	
	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline BaseSeries<PIRANHA_BASE_SERIES_TP>::~BaseSeries()
	{
        static_assert(PiranhaSeries<Derived>, "Final series class must derive from BaseSeries class.");
	}
}

#undef derived_const_cast
#undef DC

#endif
