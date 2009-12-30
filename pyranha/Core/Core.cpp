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

#include <boost/functional/hash.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/core/base_classes/named_series_def.h"
#include "../../src/core/config.h"
#include "../../src/core/mp.h"
#include "../../src/core/psym.h"
#include "../../src/core/settings.h"
#include "../args_tuple.h"
#include "../boost_python_container_conversions.h"
#include "../commons.h"
#include "../exceptions.h"
#include "../mp_classes.h"

using namespace boost::python;
using namespace piranha;
using namespace pyranha;

std::string static inline py_psym_repr(const psym &p)
{
	std::ostringstream stream;
	stream << "Symbol: '" << p.get_name() << "' - [";
	const std::size_t size = p.get_time_eval().size();
	for (std::size_t i = 0; i < size; ++i) {
		stream << p.get_time_eval()[i];
		if (i < size - 1) {
			stream << ',';
		}
	}
	stream << ']';
	return stream.str();
}

static inline std::size_t py_psym_hash(const psym &p)
{
	return boost::hash<std::string>()(p.get_name());
}

static inline void ed_set_item(eval_dict &d, const std::string &n, const double &value)
{
	d[n] = value;
}

// Instantiate the pyranha Core module.
BOOST_PYTHON_MODULE(_Core)
{
	translate_exceptions();

	// Interop between vectors of some types and Python tuples/lists.
	to_tuple_mapping<std::vector<std::string> >();
	from_python_sequence<std::vector<std::string>,variable_capacity_policy>();
	to_tuple_mapping<vector_psym>();
	from_python_sequence<vector_psym,variable_capacity_policy>();
	to_tuple_mapping<std::vector<vector_psym> >();
	to_tuple_mapping<std::vector<double> >();
	from_python_sequence<std::vector<double>,variable_capacity_policy>();
	to_tuple_mapping<std::vector<mp_rational> >();
	from_python_sequence<std::vector<mp_rational>,variable_capacity_policy>();
	to_tuple_mapping<std::vector<mp_integer> >();
	from_python_sequence<std::vector<mp_integer>,variable_capacity_policy>();

	// Expose evaluation dictionary.
	class_<eval_dict> ed("eval_dict","Evaluation dictionary.", init<>());
	ed.def("__setitem__",&ed_set_item);

	// Expose arguments tuples.
	expose_args_tuples<__PIRANHA_MAX_ECHELON_LEVEL>();

	// MP classes.
	class_<mp_rational> mpr(expose_real_mp_class<mp_rational>("rational","Multi-precision rational number."));
	mpr.def(init<const int &, const int &>());
	mpr.def(init<const mp_integer &, const mp_integer &>());
	mpr.add_property("num",&mp_rational::get_num);
	mpr.add_property("den",&mp_rational::get_den);
	mpr.def("choose", &mp_rational::choose, "Binomial coefficient (choose function).");
	class_<mp_integer> mpz(expose_real_mp_class<mp_integer>("integer","Multi-precision integer number."));
	mpz.def(init<const mp_rational &>());
	mpz.def("factorial", &mp_integer::factorial, "Factorial.");
	mpz.def("choose", &mp_integer::choose, "Binomial coefficient (choose function).");
	mpz.def("lcm", &mp_integer::lcm, "Set self to the least common multiplier of input arguments.");
	mpz.def(boost::python::self %= mp_integer());
	mpz.def(boost::python::self %= int());
	mpz.def(boost::python::self % mp_integer());
	mpz.def(boost::python::self % int());

	enum_<settings::fp_representation>("fp_representation")
	.value("scientific", settings::scientific)
	.value("decimal", settings::decimal)
	.export_values();

	class_<settings> class_setm("__settings", "Pyranha settings.", init<>());
	class_setm.add_static_property("debug", &settings::get_debug, &settings::set_debug);
	class_setm.add_static_property("used_memory", &settings::get_used_memory, "Amount of used memory in bytes.");
	class_setm.add_static_property("memory_limit", &settings::get_memory_limit, &settings::set_memory_limit);
	class_setm.add_static_property("digits", &settings::get_digits, &settings::set_digits);
	class_setm.add_static_property("fp_repr", &settings::get_fp_repr, &settings::set_fp_repr);
	class_setm.add_static_property("max_pretty_print_size", &settings::get_max_pretty_print_size, &settings::set_max_pretty_print_size);
	class_setm.add_static_property("nthread", make_function(&settings::get_nthread,return_value_policy<copy_const_reference>()),
		&settings::set_nthread);

	// Psym.
	class_<psym>("psym", "Symbol class.", init<const std::string &, const std::vector<double> &>())
		.def(init<const std::string &, const double &>())
		.def(init<const std::string &>())
		.def("__copy__", &py_copy<psym>)
		.def("__hash__", &py_psym_hash)
		.def("__repr__", &py_psym_repr)
		.def("eval", &psym::eval)
		.add_property("name", make_function(&psym::get_name,return_value_policy<copy_const_reference>()))
		.add_property("time_eval", make_function(&psym::get_time_eval,return_value_policy<copy_const_reference>()),
			&psym::set_time_eval)
		.def("list", &psym::list, "Get list of global psyms").staticmethod("list")
		.def(self == self)
		.def(self != self);
}
