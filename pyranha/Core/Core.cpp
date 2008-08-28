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
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>
#include <vector>

#include "../../src/core/base_classes/degree_truncator.h"
#include "../../src/core/base_classes/expo_truncator.h"
#include "../../src/core/base_classes/norm_truncator.h"
#include "../../src/core/config.h"
#include "../../src/core/psym.h"
#include "../../src/core/settings.h"
#include "../args_tuple.h"
#include "../cf_key_bindings.h"
#include "../commons.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace pyranha;

namespace pyranha
{
	template <class T>
	void vector_indexing(const std::string &name)
	{
		class_<std::vector<T> >((name + "_vec").c_str()).def(vector_indexing_suite<std::vector<T> >());
	}
}

// We need this wrapper because we return a reference, and this messes up the exporting of property.
static inline double load_factor_get()
{
	return settings::load_factor();
}

// Wrapper needed to emulate non-static property.
static inline psym psyms_get(const psyms &, const std::string &name)
{
	return psyms::get(name);
}

// Instantiate the pyranha Core module.
BOOST_PYTHON_MODULE(_Core)
{
	translate_exceptions();
	numerical_cfs_bindings();
	keys_bindings();
	expose_args_tuples<__PIRANHA_MAX_ECHELON_LEVEL>();

	vector_indexing<double>("double");

	// Settings.
	enum_<settings::out_format>("out_format")
	.value("plain", settings::plain)
	.value("pretty", settings::pretty)
	.value("latex", settings::latex)
	.export_values();

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
	class_setm.add_static_property("load_factor", &load_factor_get, double_set(&settings::load_factor));
	typedef void (*digits_set)(const max_fast_int &);
	class_setm.add_static_property("digits", size_t_get(&settings::digits), digits_set(&settings::digits));
	typedef settings::out_format (*format_get)();
	typedef void (*format_set)(settings::out_format);
	class_setm.add_static_property("format", format_get(&settings::format), format_set(&settings::format));
	typedef settings::fp_representation (*fp_repr_get)();
	typedef void (*fp_repr_set)(settings::fp_representation);
	class_setm.add_static_property("fp_repr", fp_repr_get(&settings::fp_repr), fp_repr_set(&settings::fp_repr));
	class_setm.add_static_property("pi_simplify", bool_get(&settings::pi_simplify), bool_set(&settings::pi_simplify));

	// Psym manager.
	class_<psyms>("__psyms", "Manager for symbols.", init<>())
	.def("__iter__", iterator<psyms, return_internal_reference<> >())
	.def("__len__", &psyms::length).staticmethod("__len__")
	.def("__repr__", &py_print_to_string<psyms>)
	.def("__getitem__", &psyms_get);

	// Psym.
	class_<psym>("psym", "Symbol class.", init<const std::string &>())
	.def(init<const std::string &, const std::string &>())
	.def(init<const std::string &, const double &>())
	.def(init<const std::string &, const double &, const double &>())
	.def(init<const std::string &, const double &, const double &, const double &>())
	.def(init<const std::string &, const double &, const double &, const double &, const double &>())
	.def(init < const std::string &, const double &, const double &, const double &, const double &,
		 const double & > ())
	.def("__copy__", &py_copy<psym>)
	.def("__repr__", &py_print_to_string<psym>)
	.def("eval", &psym::eval)
	.add_property("name", &psym::name);

	typedef void (*limit_name)(const std::string &, const max_fast_int &);
	typedef void (*limit_psym)(const piranha::psym &, const max_fast_int &);
	typedef void (*unset_void)();
	typedef void (*unset_name)(const std::string &);
	typedef void (*unset_psym)(const psym &);
	class_<base_expo_truncator>("__expo_truncator", "Exponent truncator.", init<>())
	.def("__repr__", &py_print_to_string<base_expo_truncator>)
	.def("unset", unset_void(&base_expo_truncator::unset), "Clear the list of exponent limits.")
	.def("unset", unset_name(&base_expo_truncator::unset),
		"Clear exponent limit for argument named arg1.")
	.def("unset", unset_psym(&base_expo_truncator::unset),
		"Clear exponent limit for psym arg1.").staticmethod("unset")
	.def("set", limit_name(&base_expo_truncator::set), "Set exponent limit for symbol named arg1 to integer arg2. "
		 "If arg1 does not exist, throw an error")
	.def("set", limit_psym(&base_expo_truncator::set),
		"Set exponent limit for psym arg1 to integer arg2.").staticmethod("set");

	class_<base_norm_truncator>("__norm_truncator", "Norm truncator.", init<>())
	.def("__repr__", &py_print_to_string<base_norm_truncator>)
	.def("set", &base_norm_truncator::set, "Set truncation level to 10^-arg1 of series' norm if arg1 > 0, "
		 "throw an error otherwise.").staticmethod("set")
	.def("unset", &base_norm_truncator::unset, "Disable norm-based truncation.").staticmethod("unset");

	class_<base_degree_truncator>("__degree_truncator", "Minimum degree truncator.", init<>())
	.def("__repr__", &py_print_to_string<base_degree_truncator>)
	.def("set", &base_degree_truncator::set, "Set truncation level of series minimum degree to arg1.").staticmethod("set")
	.def("unset", &base_degree_truncator::unset, "Clear minimum degree limit.").staticmethod("unset");
}
