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

#ifndef PYRANHA_ARGS_TUPLE_H
#define PYRANHA_ARGS_TUPLE_H

#include <boost/lexical_cast.hpp>
#include <boost/python/class.hpp>
#include <boost/tuple/tuple.hpp>
#include <sstream>
#include <string>

#include "../src/core/ntuple.h"
#include "../src/core/psym.h"

namespace pyranha
{
	inline void args_tuple_py_print_helper(const boost::tuples::null_type &, const std::string) {}

	template <class ArgsTuple>
	inline void args_tuple_py_print_helper(const ArgsTuple &args_tuple, std::string &out) {
		std::ostringstream stream;
		for (size_t i = 0; i < args_tuple.get_head().size(); ++i) {
			stream << "[" << i << "] " << args_tuple.get_head()[i]->name() << '\n';
		}
		out += stream.str();
		out += "-\n";
		args_tuple_py_print_helper(args_tuple.get_tail(),out);
	}

	template <int N>
	inline void expose_args_tuples() {
		boost::python::class_<typename piranha::ntuple<piranha::vector_psym_p,N>::type>
			args_tuple_inst((std::string("__args_tuple")+boost::lexical_cast<std::string>(N)+"__").c_str(),
			(std::string("Tuple of ")+boost::lexical_cast<std::string>(N)+" arguments vectors.").c_str());
		expose_args_tuples<N-1>();
	}

	template <>
	inline void expose_args_tuples<0>() {}
}

#endif
