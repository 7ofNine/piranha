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

#include "../pyranha.h"
#include "../../src/lnp.h"
#ifdef _PIRANHA_TASS
  #include "../../src/tass17/tass17.h"
#endif

BOOST_PYTHON_MODULE(_Lnp)
{
  // TASS
#ifdef _PIRANHA_TASS

  class_<tass17>("tass17","TASS theory, version 1.7.",init<>())
  .def("load", &tass17::load,"Load series into memory.")
  .staticmethod("load")
  .def("status", &tass17::status,"Print to screen TASS status.")
  .staticmethod("status")
  .def("add_delta_lambdas", &tass17::add_delta_lambdas,"Correct series with long term perturbations.")
  .staticmethod("add_delta_lambdas")
  .def("m0", &tass17::m0,return_value_policy<copy_const_reference>(),"Get Saturn's mass in Sun mass units.")
  .staticmethod("m0")
  .def("m6", &tass17::m6,return_value_policy<copy_const_reference>(),"Get inverse of Titan's mass in Saturn mass units.")
  .staticmethod("m6")
  .def("lambda4", &tass17::lambda4,return_value_policy<copy_const_reference>(),"Get Dione's lambda.")
  .staticmethod("lambda4")
  .def("lambda6", &tass17::lambda6,return_value_policy<copy_const_reference>(),"Get Titan's lambda.")
  .staticmethod("lambda6")
  .def("r6", &tass17::r6,"Calculate Titan's radius.")
  .staticmethod("r6")
  .def("p6", &tass17::p6,return_value_policy<copy_const_reference>(),"Get Titan's p.")
  .staticmethod("p6")
  .def("z6", &tass17::z6,return_value_policy<copy_const_reference>(),"Get Titan's z.")
  .staticmethod("z6")
  .def("zeta6", &tass17::zeta6,return_value_policy<copy_const_reference>(),"Get Titan's (greek) zeta.")
  .staticmethod("zeta6")
  .def("e", &tass17::e,"Calculate eccentricity from elliptic orbital element z.")
  .staticmethod("e")
  .def("a", &tass17::a,"Calculate semi-major axis a from elliptic orbital element p.")
  .staticmethod("a")
  .def("eiM", &tass17::eiM,"Calculate complex exponential of mean mean motion M.")
  .staticmethod("eiM")
  .def("vienne_r", &tass17::vienne_r,"Calculate radius using Vienne's FORTRAN routine.")
  .staticmethod("vienne_r")
  ;

  class_ <tc_vienne_r6>
  ((std::string("tc_vienne_r6")).c_str(),
   init<tc_vienne_r6::b_type,double,double,size_t>())
  .def("sigma",&tc_vienne_r6::sigma,return_value_policy<copy_const_reference>())
  .def("max_error",&tc_vienne_r6::max_error,
       return_value_policy<copy_const_reference>())
  .def("time",&tc_vienne_r6::time,return_value_policy<copy_const_reference>())
  .def("hs",&tc_vienne_r6::hs,return_value_policy<copy_const_reference>())
  .def("hs_computed",&tc_vienne_r6::hs_computed,
       return_value_policy<copy_const_reference>())
  .def("error",&tc_vienne_r6::error)
  .def("gnuplot_save",&tc_vienne_r6::gnuplot_save)
  .def("size",&tc_vienne_r6::size)
  ;
#endif

  class_<lnp> inst=ps_basic_instantiation<lnp>("lnp","Numerical Poisson series class, "
                   "trigonometric lists version.");
  ps_instantiate_real_specifics(inst);
}
