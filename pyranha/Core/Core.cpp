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

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <string>
#include <vector>

#include "../../src/core/base_classes/degree_truncator.h"
#include "../../src/core/base_classes/norm_truncator.h"
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

// Instantiate the pyranha Core module.
BOOST_PYTHON_MODULE(_Core)
{
	translate_exceptions();

	// Interop between vectors of some types and Python tuples/lists.
	to_tuple_mapping<vector_psym>();
	from_python_sequence<vector_psym,variable_capacity_policy>();
	to_tuple_mapping<std::vector<vector_psym> >();
	to_tuple_mapping<std::vector<double> >();
	from_python_sequence<std::vector<double>,variable_capacity_policy>();
	to_tuple_mapping<std::vector<mp_rational> >();
	from_python_sequence<std::vector<mp_rational>,variable_capacity_policy>();
	to_tuple_mapping<std::vector<mp_integer> >();
	from_python_sequence<std::vector<mp_integer>,variable_capacity_policy>();

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
	mpz.def(boost::python::self %= mp_integer());
	mpz.def(boost::python::self %= int());
	mpz.def(boost::python::self % mp_integer());
	mpz.def(boost::python::self % int());

	enum_<settings::fp_representation>("fp_representation")
	.value("scientific", settings::scientific)
	.value("decimal", settings::decimal)
	.export_values();

	typedef bool (*bool_get)();
	typedef void (*bool_set)(const bool &);
	typedef size_t (*size_t_get)();
	typedef void (*size_t_set)(const size_t &);
	typedef double (*double_get)();
	typedef void (*double_set)(const double &);
	class_<settings> class_setm("__settings", "Pyranha settings.", init<>());
	class_setm.add_static_property("debug", bool_get(&settings::debug), bool_set(&settings::debug));
	class_setm.add_static_property("used_memory", &settings::used_memory, "Amount of used memory in bytes.");
	class_setm.add_static_property("memory_limit", size_t_get(&settings::memory_limit),
		size_t_set(&settings::memory_limit));
	class_setm.add_static_property("load_factor", make_function(&settings::get_load_factor,return_value_policy<copy_const_reference>()),
		&settings::set_load_factor);
	typedef void (*digits_set)(const int &);
	class_setm.add_static_property("digits", size_t_get(&settings::digits), digits_set(&settings::digits));
	typedef settings::fp_representation (*fp_repr_get)();
	typedef void (*fp_repr_set)(settings::fp_representation);
	class_setm.add_static_property("fp_repr", fp_repr_get(&settings::fp_repr), fp_repr_set(&settings::fp_repr));
	class_setm.add_static_property("max_pretty_print_size", &settings::get_max_pretty_print_size, &settings::set_max_pretty_print_size);

	// Psym.
	class_<psym>("psym", "Symbol class.", init<const std::string &, const std::vector<double> &>())
		.def(init<const std::string &, const double &>())
		.def(init<const std::string &>())
		.def("__copy__", &py_copy<psym>)
		.def("__repr__", &py_print_to_string<psym>)
		.def("eval", &psym::eval)
		.add_property("name", make_function(&psym::get_name,return_value_policy<copy_const_reference>()))
		.add_property("time_eval", make_function(&psym::get_time_eval,return_value_policy<copy_const_reference>()),
			&psym::set_time_eval)
		.def("list", &psym::list, "Get list of global psyms").staticmethod("list");

	class_<norm_truncator>("__norm_truncator", "Norm truncator.", init<>())
	.def("__repr__", &py_print_to_string<norm_truncator>)
	.def("set", &norm_truncator::set, "Set truncation level to 10^-arg1 of series' norm if arg1 > 0, "
		 "throw an error otherwise.").staticmethod("set")
	.def("unset", &norm_truncator::unset, "Disable norm-based truncation.").staticmethod("unset");

	typedef void (*deg_set)(const int &);
	typedef void (*p_deg_set)(const vector_psym &, const int &);
	class_<degree_truncator>("__degree_truncator", "Minimum degree truncator.", init<>())
	.def("__repr__", &py_print_to_string<degree_truncator>)
	.def("set", deg_set(&degree_truncator::set), "Set truncation level of series' minimum degree to arg1.")
	.def("set", p_deg_set(&degree_truncator::set), "Set truncation level of series' partial minimum degree to arg2, "
		"relatively to list of psyms psym1.").staticmethod("set")
	.def("unset", &degree_truncator::unset, "Clear minimum degree limit.").staticmethod("unset");
}
