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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>
#include <vector>

#include "../../src/core/base_classes/expo_truncator.h"
#include "../../src/core/base_classes/norm_truncator.h"
#include "../../src/core/psym.h"
#include "../../src/core/settings.h"
#include "../../src/core/stats.h"
#include "../../src/core/stream_manager.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace pyranha;

namespace pyranha
{
  template <class T>
    void vector_indexing(const std::string &name)
  {
    class_<std::vector<T> >((name+"_vec").c_str()).def(vector_indexing_suite<std::vector<T> >());
  }
}

// Instantiate the pyranha Core module.
BOOST_PYTHON_MODULE(_Core)
{
  translate_exceptions();
  vector_indexing<double>("double");

  // Settings.
  typedef const bool &(*debug_get)();
  typedef void (*debug_set)(const bool &);
  class_<settings> class_setm("settings","Settings for Pyranha.",init<>());
  class_setm.def("debug",debug_get(&settings::debug),return_value_policy<copy_const_reference>(),
    "Get value of the debug flag.");
  class_setm.def("debug",debug_set(&settings::debug),"Set value of the debug flag.").staticmethod("debug");
//   class_setm.def("debug",debug_get(&settings_manager::debug),return_value_policy<copy_const_reference>(),
//     "Get value of the debug flag").staticmethod("debug");
//   class_setm.def("load_factor", &settings_manager::load_factor,return_value_policy<copy_const_reference>(),
//     "Get value of maximum load factor for hashed containers.");
//   class_setm.def("numerical_zero", &settings_manager::numerical_zero,return_value_policy<copy_const_reference>(),
//     "Get value of numerical zero.");
//   class_setm.def("theories_path", &settings_manager::theories_path,return_value_policy<copy_const_reference>(),
//     "Get search path for theories of motion's data files.");
//   class_setm.def("default_theories_path", &settings_manager::default_theories_path,return_value_policy<copy_const_reference>(),
//     "Get default search path for theories of motion's data files.");
//   class_setm.def("version", &settings_manager::version,
//     "Get Piranha's version.", return_value_policy<copy_const_reference>());
//   class_setm.def("display_progress", &settings_manager::display_progress,
//     "Display progress bar?");
//   class_setm.def("set_display_progress", &settings_manager::set_display_progress,
//     "Set to true to enable the display of progress bar.");
//   class_setm.def("mp_default_prec", &settings_manager::mp_default_prec,
//     "Get least default precision of mp floating point in bits.");
//   class_setm.def("set_mp_default_prec", &settings_manager::set_mp_default_prec,
//     "Set least default precision of mp floating point in bits.");
//   class_setm.def("set_load_factor", &settings_manager::set_load_factor,
//     "Set value of maximum load factor for hashed containers.");
//   class_setm.def("set_theories_path", &settings_manager::set_theories_path,
//     "Set search path for theories of motion's data files.");
//   class_setm.staticmethod("load_factor");
//   class_setm.staticmethod("set_load_factor");
//   class_setm.staticmethod("numerical_zero");
//   class_setm.staticmethod("theories_path");
//   class_setm.staticmethod("default_theories_path");
//   class_setm.staticmethod("set_theories_path");
//   class_setm.staticmethod("version");
//   class_setm.staticmethod("display_progress");
//   class_setm.staticmethod("set_display_progress");
//   class_setm.staticmethod("mp_default_prec");
//   class_setm.staticmethod("set_mp_default_prec");

  // Stream manager.
  enum_<stream_manager::out_format>("out_format")
    .value("plain", stream_manager::plain)
    .value("latex", stream_manager::latex)
    .export_values();

  enum_<stream_manager::fp_representation>("fp_representation")
    .value("scientific", stream_manager::scientific)
    .value("decimal", stream_manager::decimal)
    .export_values();

  class_<stream_manager> class_sm("stream_manager","Set up stream output.",no_init);
  class_sm.def("digits",&stream_manager::digits,"Get number of digits used in output.");
  class_sm.def("min_digits",&stream_manager::min_digits,"Get minimum number of digits used in output.");
  class_sm.def("max_digits",&stream_manager::max_digits,"Get maximum number of digits used in output.");
  class_sm.def("set_digits",&stream_manager::set_digits,"Set number of digits used in output.");
  class_sm.def("format",&stream_manager::format,"Get stream output format.");
  class_sm.def("set_format",&stream_manager::set_format,"Set stream output format.");
  class_sm.def("set_fp_rep", &stream_manager::set_fp_rep,
    "Set floating-point representation.");
  class_sm.def("fp_rep",&stream_manager::fp_rep,"Get floating-point representation.");
  class_sm.staticmethod("digits");
  class_sm.staticmethod("min_digits");
  class_sm.staticmethod("max_digits");
  class_sm.staticmethod("set_digits");
  class_sm.staticmethod("format");
  class_sm.staticmethod("set_format");
  class_sm.staticmethod("fp_rep");
  class_sm.staticmethod("set_fp_rep");

  // Stats.
  class_<stats>("stats","Piranha-specific statistics.",no_init)
    .def("pack_ratio",&stats::pack_ratio)
    .staticmethod("pack_ratio");

  // Psym manager.
  class_<psym_manager>("__psym_manager","Manager for symbols.",init<>())
    .def("__iter__",iterator<psym_manager,return_internal_reference<> >()).staticmethod("__iter__")
    .def("__len__",&psym_manager::length).staticmethod("__len__")
    .def("__repr__",&psym_manager::print_to_string).staticmethod("__repr__");

  // Psym.
  class_<psym>("psym","Symbol class.",init<const std::string &>())
    .def(init<const std::string &, const std::string &>())
    .def(init<const std::string &, const double &>())
    .def(init<const std::string &, const double &, const double &>())
    .def(init<const std::string &, const double &, const double &, const double &>())
    .def(init<const std::string &, const double &, const double &, const double &, const double &>())
    .def(init<const std::string &, const double &, const double &, const double &, const double &,
    const double &>())
    .def("__copy__",&psym::copy)
    .def("__repr__",&psym::print_to_string);

  typedef void (*limit_name)(const std::string &, const max_fast_int &);
  typedef void (*limit_psym)(const piranha::psym &, const max_fast_int &);
  class_<base_expo_truncator>("__expo_truncator","Exponent truncator.",init<>())
    .def("__repr__",&base_expo_truncator::print_to_string).staticmethod("__repr__")
    .def("clear_all",&base_expo_truncator::clear_all,"Clear list of exponent limits.").staticmethod("clear_all")
    .def("clear",&base_expo_truncator::clear,"Clear exponent limit for argument named arg1.").staticmethod("clear")
    .def("limit",limit_name(&base_expo_truncator::limit),"Set exponent limit for symbol named arg1 to integer arg2. "
      "If arg1 does not exist, throw an error")
    .def("limit",limit_psym(&base_expo_truncator::limit),"Set exponent limit for psym arg1 to integer arg2.").staticmethod("limit");

  class_<base_norm_truncator>("__norm_truncator","Norm truncator.",init<>())
    .def("__repr__",&base_norm_truncator::print_to_string).staticmethod("__repr__")
    .def("set",&base_norm_truncator::set,"Set truncation level of series norm to 10^-arg1 if arg1 > 0, to 0 if arg1 == 0 "
      "and throw an error otherwise.").staticmethod("set");
}
