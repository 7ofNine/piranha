/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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
#include <boost/python/return_value_policy.hpp>

#include "../../src/astro.h"
#include "../../src/phase_list.h"
#include "../../src/psymbol.h"
#include "../../src/stats.h"
#include "../../src/stream_manager.h"

using namespace boost::python;
using namespace piranha;

// Instantiate the pyranha Core module.
BOOST_PYTHON_MODULE(_Core)
{
  // Astronomical class instantiation.
  class_<astro> class_astro("astro","Useful astronomical functions and constants.",no_init);
  class_astro.def("G",&astro::G,return_value_policy<copy_const_reference>(),
                  "Universal gravitational constant, from Standish (1995).");
  class_astro.def("k",&astro::k,return_value_policy<copy_const_reference>(),
                  "Gaussian gravitational constant.");
  class_astro.def("eps_0",&astro::eps_0,return_value_policy<copy_const_reference>(),
                  "Obliquity of ecliptic at J2000, from Standish (1995).");
  class_astro.def("J2000dot0",&astro::J2000dot0,return_value_policy<copy_const_reference>(),
                  "J2000.0 epoch in Julian days.");
  class_astro.def("J1980dot0",&astro::J1980dot0,return_value_policy<copy_const_reference>(),
                  "J1980.0 epoch in Julian days.");
  class_astro.def("JD_per_JY",&astro::JD_per_JY,return_value_policy<copy_const_reference>(),
                  "Julian Days per Julian Year.");
  class_astro.def("seconds_per_JY",&astro::seconds_per_JY,return_value_policy<copy_const_reference>(),
                  "Seconds per Julian Year.");
  class_astro.def("AU",&astro::AU,return_value_policy<copy_const_reference>(),
                  "Astronomical Unit, from http://ssd.jpl.nasa.gov/?astro.");
  // Static methods instantiations.
  class_astro.staticmethod("G");
  class_astro.staticmethod("k");
  class_astro.staticmethod("eps_0");
  class_astro.staticmethod("J2000dot0");
  class_astro.staticmethod("J1980dot0");
  class_astro.staticmethod("JD_per_JY");
  class_astro.staticmethod("seconds_per_JY");
  class_astro.staticmethod("AU");
  // Astronomical functions.
  class_astro.def("JD_to_elp2000",&astro::JD_to_elp2000);
  class_astro.def("kep_cosE",&astro::kep_cosE<double>,"Solve Kepler's equation for cosE.");
  class_astro.def("sph_to_x",&astro::sph_to_x,"Convert spherical coordinates into x coordinate.");
  class_astro.def("sph_to_y",&astro::sph_to_y,"Convert spherical coordinates into y coordinate.");
  class_astro.def("sph_to_z",&astro::sph_to_z,"Convert spherical coordinates into z coordinate.");
  class_astro.staticmethod("JD_to_elp2000");
  class_astro.staticmethod("kep_cosE");
  class_astro.staticmethod("sph_to_x");
  class_astro.staticmethod("sph_to_y");
  class_astro.staticmethod("sph_to_z");


  // Instantiate mathematical functions.
  class_<math> class_math("math","Pyranha mathematical functions for double precision numbers.",no_init);
  class_math.def("norm",math::norm<double>,"Norm.");
  class_math.def("natural_pow",math::natural_pow<double>,"Natural power.");
  class_math.def("Pnm",math::Pnm<double>,"Legendre function of the first kind - Pnm(cos(theta)).");
  class_math.def("complexp",math::complexp<double>,"Complex exponential.");
  class_math.def("cosine",math::cosine<double>,"Cosine.");
  class_math.def("sine",math::sine<double>,"Sine.");
  class_math.def("Ynm",math::Ynm<double>,"Non-normalized spherical harmonic.");
  class_math.def("wig_rot",math::wig_rot<double>,"Wigner rotation theorem for spherical harmonics.");
  class_math.staticmethod("norm");
  class_math.staticmethod("natural_pow");
  class_math.staticmethod("cosine");
  class_math.staticmethod("sine");
  class_math.staticmethod("complexp");
  class_math.staticmethod("Pnm");
  class_math.staticmethod("Ynm");
  class_math.staticmethod("wig_rot");


  // Settings.
  class_<settings_manager> class_setm("settings_manager","Manage piranha-specific settings.",no_init);
  class_setm.def("set_prec", &settings_manager::set_prec,
                 "Set precision of mathematical operations.");
  class_setm.def("prec", &settings_manager::prec,return_value_policy<copy_const_reference>(),
                 "Get precision of mathematical operations.");
  class_setm.def("numerical_zero", &settings_manager::numerical_zero,return_value_policy<copy_const_reference>(),
                 "Get value of numerical zero.");
  class_setm.def("theories_path", &settings_manager::theories_path,return_value_policy<copy_const_reference>(),
                 "Get search path for theories of motion's data files.");
  class_setm.def("default_theories_path", &settings_manager::default_theories_path,return_value_policy<copy_const_reference>(),
                 "Get default search path for theories of motion's data files.");
  class_setm.def("set_theories_path", &settings_manager::set_theories_path,
                 "Set search path for theories of motion's data files.");
  class_setm.staticmethod("set_prec");
  class_setm.staticmethod("prec");
  class_setm.staticmethod("numerical_zero");
  class_setm.staticmethod("theories_path");
  class_setm.staticmethod("default_theories_path");
  class_setm.staticmethod("set_theories_path");


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
  class_sm.def("data_separator",&stream_manager::data_separator,
               return_value_policy<copy_const_reference>(),"Get string used as data separator.");
  class_sm.def("format",&stream_manager::format,"Get stream output format.");
  class_sm.def("set_format",&stream_manager::set_format,"Set stream output format.");
  class_sm.def("set_fp_rep", &stream_manager::set_fp_rep,
               "Set floating-point representation.");
  class_sm.def("fp_rep",&stream_manager::fp_rep,"Get floating-point representation.");
  class_sm.staticmethod("digits");
  class_sm.staticmethod("min_digits");
  class_sm.staticmethod("max_digits");
  class_sm.staticmethod("set_digits");
  class_sm.staticmethod("data_separator");
  class_sm.staticmethod("format");
  class_sm.staticmethod("set_format");
  class_sm.staticmethod("fp_rep");
  class_sm.staticmethod("set_fp_rep");


  // Stats.
  class_<stats>("stats","Piranha-specific statistics.",no_init)
  .def("pack_ratio",&stats::pack_ratio)
  .staticmethod("pack_ratio");


  // Symbols.
  class_<psymbol_manager>("psymbol_manager","Manager for psymbols.",no_init)
  .def("put", &psymbol_manager::put,"Show registered symbols.")
  .staticmethod("put")
  ;


  // List of phases.
  // FIXME: expose method to manipulate them?
  class_<phase_list>("phase_list","List of phases.",init<std::string>());

  // Psymbols.
  class_<psymbol>("psymbol","Symbol class.")
  .def(init<const std::string &, const double &>())
  .def(init<const std::string &, const double &, const double &>())
  .def(init<const std::string &, const double &, const double &, const double &>())
  .def(init<const std::string &, const double &, const double &, const double &, const double &>())
  .def(init<const std::string &, const double &, const double &, const double &, const double &,
     const double &>())
  .def("__copy__",&psymbol::copy)
  .def("put",&psymbol::put)
  .def("name",&psymbol::name,return_value_policy<copy_const_reference>())
  .def("phase",&psymbol::phase)
  .def("freq",&psymbol::freq);
}
