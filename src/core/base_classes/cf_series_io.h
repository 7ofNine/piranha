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

#ifndef PIRANHA_CF_SERIES_IO_H
#define PIRANHA_CF_SERIES_IO_H

#include <boost/algorithm/string.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace piranha
{
	/// Constructor from string.
	/**
	  * The whole series is stored into a string when using it as coefficient in another series.
	  */
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void CfSeries<__PIRANHA_CF_SERIES_TP>::constructFromString(const std::string &input, const ArgsTuple &argsTuple)
	{
		typedef typename Derived::TermType TermType;
		const char separator = Derived::separator;

		std::string str(input);
		// Remove extra spaces.
		boost::trim(str);
		// Split into single terms.
		std::vector<std::string> terms;
		boost::split(terms, str, boost::is_any_of(std::string(1, separator)));

		const std::size_t length = terms.size();
		for (std::size_t i = 0; i < length; ++i) 
        {
			// Build the term from the string.
			TermType term(terms[i], argsTuple);
			// Insert it.
			derived_cast->insert(term, argsTuple);
		}
	}


	/// Print in plain mode.
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void CfSeries<__PIRANHA_CF_SERIES_TP>::printPlain(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		derived_const_cast->printTermsPlain(stream, argsTuple);
	}


	/// Print in pretty mode.
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void CfSeries<__PIRANHA_CF_SERIES_TP>::printPretty(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		if (derived_const_cast->length() > 1) 
        {
			stream << '(';
		}

		derived_const_cast->printTermsPretty(stream, argsTuple);
		
        if (derived_const_cast->length() > 1) 
        {
			stream << ')';
		}
	}


	/// Print in tex mode.
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void CfSeries<__PIRANHA_CF_SERIES_TP>::printTex(std::ostream &stream, const ArgsTuple &argsTuple) const
	{
		if (derived_const_cast->length() > 1) 
        {
			stream << "\\left(";
		}

		derived_const_cast->printTermsTEX(stream, argsTuple);
		
        if (derived_const_cast->length() > 1) 
        {
			stream << "\\right)";
		}
	}
}

#endif
