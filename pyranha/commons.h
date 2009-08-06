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

#ifndef PYRANHA_COMMONS_H
#define PYRANHA_COMMONS_H

#include <boost/tuple/tuple.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "../src/core/psym.h"

namespace pyranha
{
	template <class T>
	inline T py_copy(const T &x) {
		return T(x);
	}

	template <class T>
	std::string py_print_to_string(const T &origin)
	{
		std::ostringstream stream;
		origin.print(stream);
		return stream.str();
	}

	template <class T>
	std::string py_print_to_string_tex(const T &origin)
	{
		std::ostringstream stream;
		origin.print_tex(stream);
		return stream.str();
	}

	template <class T>
	std::string py_print_to_string_plain(const T &origin)
	{
		std::ostringstream stream;
		origin.print_plain(stream);
		return stream.str();
	}

	template <class Tuple>
	struct py_tuple_to_vector
	{
		template <class Vector>
		static void run(Vector &retval, const Tuple &tuple)
		{
			retval.push_back(tuple.get_head());
			py_tuple_to_vector<typename Tuple::tail_type>::template run(retval,tuple.get_tail());
		}
	};

	template <>
	struct py_tuple_to_vector<boost::tuples::null_type>
	{
		template <class Vector>
		static void run(Vector &, const boost::tuples::null_type &) {}
	};

	template <class Series>
	inline std::vector<piranha::vector_psym> py_series_arguments(const Series &s)
	{
		std::vector<piranha::vector_psym> retval;
		py_tuple_to_vector<typename Series::args_tuple_type>::run(retval,s.arguments());
		return retval;
	}
}

#endif
