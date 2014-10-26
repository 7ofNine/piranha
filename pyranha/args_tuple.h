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

#ifndef PYRANHA_argsTuple_H
#define PYRANHA_argsTuple_H

#include <boost/lexical_cast.hpp>
#include <boost/python/class.hpp>
#include <boost/tuple/tuple.hpp>
#include <sstream>
#include <string>

#include "../src/core/ntuple.h"
#include "../src/core/Psym.h"

namespace pyranha
{
	inline void argsTuple_py_print_helper(const boost::tuples::null_type &, const std::string) {}

	template <class ArgsTuple>
	inline void argsTuple_py_print_helper(const ArgsTuple &argsTuple, std::string &out)
	{
		std::ostringstream stream;
		for (size_t i = 0; i < argsTuple.get_head().size(); ++i) {
			stream << i << ' ' << argsTuple.get_head()[i].get_name() << '\n';
		}
		out += stream.str();
		argsTuple_py_print_helper(argsTuple.get_tail(), out);
	}

	template <class ArgsTuple>
	inline std::string py_argsTuple_repr(const ArgsTuple &argsTuple)
	{
		std::string retval;
		argsTuple_py_print_helper(argsTuple, retval);
		return retval;
	}

	template <int N>
	inline void expose_argsTuples()
	{
		typedef typename piranha::NTuple<piranha::VectorPsym, N>::Type ArgsTupleType;
		boost::python::class_<ArgsTupleType>
		argsTuple_inst((std::string("__base_argsTuple") + boost::lexical_cast<std::string>(N) + "__").c_str(),
						(std::string("Tuple of ") + boost::lexical_cast<std::string>(N) + " arguments vectors.").c_str());
		argsTuple_inst.def(boost::python::init<const ArgsTupleType &>());
		argsTuple_inst.def("__repr__", &py_argsTuple_repr<ArgsTupleType>);
		expose_argsTuples < N - 1 > ();
	}

	template <>
	inline void expose_argsTuples<0>() {}
}

#endif
