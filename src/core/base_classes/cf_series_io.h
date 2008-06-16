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
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::construct_from_string(const std::string &str_, const ArgsTuple &args_tuple)
	{
		typedef typename Derived::term_type term_type;
		const char separator = Derived::separator;
		std::string str(str_);
		// Remove extra spaces.
		boost::trim(str);
		// Split into single terms.
		std::vector<std::string> vs;
		boost::split(vs, str, boost::is_any_of(std::string(1, separator)));
		const size_t length = vs.size();
		for (size_t i = 0; i < length; ++i) {
			try {
				// Try to build the term from the string.
				// Here we check that the term is insertable, i.e., the number of arguments is compatible
				// with the provided args tuple.
				term_type term(vs[i], args_tuple);
				// TODO: most likely this throw can be moved inside the main insert function of base_series,
				// as soon as we make sure that it won't impact performance too much.
				if (!term.is_insertable(args_tuple)) {
					throw term_not_insertable("Term not insertable in cf_series.");
				}
				derived_cast->insert(term, args_tuple);
			} catch (const bad_input &bi) {
				std::cout << bi.what() << std::endl;
			} catch (const term_not_insertable &tni) {
				std::cout << tni.what() << std::endl;
			}
		}
	}

	/// Print in plain mode.
	template <__PIRANHA_CF_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void cf_series<__PIRANHA_CF_SERIES_TP>::print_plain(std::ostream &stream, const ArgsTuple &args_tuple) const
	{
		derived_const_cast->print_terms_plain(stream, args_tuple, -1);
	}
}

#endif
