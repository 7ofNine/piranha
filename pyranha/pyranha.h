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

#ifndef PYRANHA_PYRANHA_H
#define PYRANHA_PYRANHA_H

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/self.hpp>
#include <complex>
#include <exception>
#include <vector>

#include "../src/piranha.h"

using namespace boost::python;
using namespace piranha;

#include "stl_containers.h"
#include "exceptions.h"

/// Instantiation of spectral comparison.
// template <class T>
//   class_<sc<T> > sc_instatiation()
// {
//   class_<sc<T> > retval("sc",init<T,T,bool>());
//   retval.def(init<T,T>());
//   retval.def("size",&sc<T>::size);
//   retval.def("is_relative",&sc<T>::is_relative);
//   retval.def("diffs",&sc<T>::diffs,return_value_policy<copy_const_reference>());
//   retval.def("gnuplot_save", &sc<T>::gnuplot_save);
//   return retval;
// }

// template <class T, bool parallel>
//   void range_evaluator_instantiation(const std::string &name, const std::string &descr)
// {
//   class_<range_evaluator<T,parallel> > range_eval(name.c_str(),descr.c_str(),
//       init<const T &,const double &, const double &, const int &>());
//   range_eval.def("values", &range_evaluator<T,parallel>::values, return_internal_reference<1>());
//   range_eval.def("times", &range_evaluator<T,parallel>::times, return_internal_reference<1>());
// }

// template <class T>
//   void range_evaluator_instantiation(const std::string &name, const std::string &descr)
// {
//   range_evaluator_instantiation<T,compile_switches::use_tbb>(name,descr);
// }

/// Instantiation of common functions for time comparison.
// Here T is the tc of operation<series_type>.
// template <class T>
//   void tc_common_instantiation(const std::string &name, class_<T> &time_c)
// {
//   time_c.def("sigma",&T::sigma,return_value_policy<copy_const_reference>());
//   time_c.def("max_error",&T::max_error,return_value_policy<copy_const_reference>());
//   time_c.def("gnuplot_save",&T::gnuplot_save);
//   range_evaluator_instantiation<T>(std::string("__range_evaluator_")+name,
//     std::string("Range evaluator for ")+name+".");
//   time_c.def("plain_eval",&T::plain_eval,return_internal_reference<1>());
//   time_c.def("manip_eval",&T::manip_eval,return_internal_reference<1>());
// }

/// Helper macro that calls tc_common_instantiation using just one parameter instead of two.
#define __tc_common_instantiation(name) tc_common_instantiation(#name,name##_inst)

/// Template for the instantiation of a ps. It will expose methods common to real and complex
/// ps, the specializations take place below.
template <class T>
  class_<T> series_basic_instantiation(const std::string &name, const std::string &description)
{
// This is a trick to help resolve overloaded methods inside classes.
  //typedef void (T::*crop_it)(const typename T::it_s_index &);
  //typedef void (T::*crop_real)(const double &);
//   typedef typename eval_type<T>::type (T::*eval_single) (const double &) const;
//     const size_t &) const;
  class_<T> inst(name.c_str(),description.c_str());
  inst.def(init<const T &>());
  inst.def(init<const std::string &>());
  inst.def(init<const int &>());
  inst.def(init<const double &>());
  inst.def("__copy__",&T::copy);
  inst.def("__repr__",&T::print_to_string);
//   inst.def("__iter__", iterator<T,return_internal_reference<> >());
  inst.def("__len__",&T::length);
//   inst.def("begin", &T::begin);
//   inst.def("end", &T::end);
//   inst.def("address", &T::address);
  inst.def("save_to",&T::save_to, "Save series to file.");
//   inst.def("put_phases_freqs", put_phases_freqs_noargs(&T::put_phases_freqs));
//   inst.def("put_phases_freqs", put_phases_freqs_n(&T::put_phases_freqs));
  inst.def("length",&T::length);
  inst.def("eval",&T::eval);
//   inst.def("cf_width", &T::cf_width);
//   inst.def("trig_width", &T::trig_width);
//   inst.def("is_empty", &T::is_empty);
//   inst.def("g_norm", &T::g_norm);
//   inst.def("footprint", &T::footprint);
//   inst.def("checkup", &T::checkup);
  //inst.def("crop", crop_real(&T::crop));
//   inst.def("t_eval", t_eval_single(&T::t_eval));
//   inst.def("t_eval_brute", &T::t_eval_brute);
//   inst.def("mean", mean_def(&T::mean));
//   inst.def("mean", mean_n(&T::mean));
  inst.def("swap",&T::swap);
// NOTICE: the order seems important here, if we place *=int before *=double we
// will get just *=double in Python. Go figure...
// Assignments.
//   inst.def(self=int());
//   inst.def(self=double());
//   inst.def(self=self);
// Addition and subtraction.
  inst.def(self+=int());
  inst.def(self+=double());
  inst.def(self+=self);
  inst.def(self+int());
  inst.def(int()+self);
  inst.def(self+double());
  inst.def(double()+self);
  inst.def(self+self);
  inst.def(self-=int());
  inst.def(self-=double());
  inst.def(self-=self);
  inst.def(self-int());
  inst.def(int()-self);
  inst.def(self-double());
  inst.def(double()-self);
  inst.def(self-self);
// Multiplication.
  inst.def(self*=int());
  inst.def(self*=double());
  inst.def(self*=self);
  inst.def(self*int());
  inst.def(int()*self);
  inst.def(self*double());
  inst.def(double()*self);
  inst.def(self*self);
// Division.
  inst.def(self/=int());
  inst.def(self/=double());
  inst.def(self/int());
  inst.def(self/double());
// Instantiate spectral comparison.
  //sc_instatiation<T>();
// Instantiate common time comparisons.
//   class_<tc_equal<T> > tc_equal_inst("tc_equal",
//     init<typename tc_equal<T>::b_type,double,double,size_t,T>());
//   __tc_common_instantiation(tc_equal);
//   class_<tc_mult<T> > tc_mult_inst("tc_mult",
//     init<typename tc_mult<T>::b_type,double,double,size_t,T,T>());
//   __tc_common_instantiation(tc_mult);
//   class_<tc_insert_phases<T> > tc_insert_phases_inst("tc_insert_phases",
//     init<typename tc_insert_phases<T>::b_type,double,double,size_t,phase_list,T>());
//   __tc_common_instantiation(tc_insert_phases);
// Range evaluator for series.
//   range_evaluator_instantiation<T>("range_evaluator",std::string("Evaluate ")+name+
//     " over a time range.");
  return inst;
}

template <class T>
  void series_trigonometric_instantiation(class_<T> &inst)
{
  inst.def("cos",&T::cos);
  inst.def("sin",&T::sin);
}

template <class T>
  void series_pow_instantiation(class_<T> &inst)
{
  typedef T (T::*pow_unary)(const double &) const;
  inst.def("__pow__",pow_unary(&T::pow));
}

template <class T>
  void series_psymbol_instantiation(class_<T> &inst)
{
  inst.def(init<const psymbol &>());
}

template <class T>
  void ps_instantiate_differential_specifics(class_<T> &inst)
{
  inst.def("partial", &T::partial);
}

template <class T>
  void fourier_specifics(class_<T> &inst)
{

}

// template <class T>
//   void ps_instantiate_real_specifics(class_<T> &real)
// {
//   typedef T real_ps;
//   real.def(init<const psymbol &, psymbol::type>());
//   real.def("complexp", &real_ps::complexp);
//   real.def("cosine", &real_ps::cosine);
//   real.def("sine", &real_ps::sine);
//   real.def("pow", &real_ps::pow);
// // External functions.
// //   def("kep_cosE",&astro::kep_cosE<real_ps>,"Solve Kepler's equation for cosE.");
// // // NOTE: which functions does it make sense to keep here?
// //   def("Pnm",&math::Pnm<real_ps>,"Legendre function of the first kind - Pnm(cos(theta)).");
// //   def("Ynm",&math::Ynm<real_ps>,"Non-normalized spherical harmonic.");
// //   def("wig_rot",&math::wig_rot<real_ps>,"Wigner rotation theorem for spherical harmonics.");
//   class_<tc_complexp<real_ps> > tc_complexp_inst("tc_complexp",
//     init<typename tc_complexp<real_ps>::b_type,double,double,size_t,real_ps>());
//   __tc_common_instantiation(tc_complexp);
//   class_<tc_cosine<real_ps> > tc_cosine_inst("tc_cosine",
//     init<typename tc_cosine<real_ps>::b_type,double,double,size_t,real_ps>());
//   __tc_common_instantiation(tc_cosine);
//   class_<tc_sine<real_ps> > tc_sine_inst("tc_sine",
//     init<typename tc_sine<real_ps>::b_type,double,double,size_t,real_ps>());
//   __tc_common_instantiation(tc_sine);
//   class_<tc_Pnm<real_ps> > tc_Pnm_inst("tc_Pnm",
//     init<typename tc_Pnm<real_ps>::b_type,double,double,size_t,int,int,real_ps>());
//   __tc_common_instantiation(tc_Pnm);
//   class_<tc_Ynm<real_ps> > tc_Ynm_inst("tc_Ynm",
//     init<typename tc_Ynm<real_ps>::b_type,double,double,size_t,int,int,real_ps,real_ps>());
//   __tc_common_instantiation(tc_Ynm);
//   class_<tc_wig_rot<real_ps> > tc_wig_rot_inst("tc_wig_rot",
//     init<typename tc_wig_rot<real_ps>::b_type,double,double,size_t,int,int,real_ps,real_ps,
//     real_ps,real_ps,real_ps>());
//   __tc_common_instantiation(tc_wig_rot);
//   class_<tc_pow<real_ps> > tc_pow_inst("tc_pow",
//     init<typename tc_pow<real_ps>::b_type,double,double,size_t,double,real_ps>());
//   __tc_common_instantiation(tc_pow);
//   class_<tc_add_ps_to_arg<real_ps> > tc_add_ps_to_arg_inst("tc_add_ps_to_arg",
//     init<typename tc_add_ps_to_arg<real_ps>::b_type,double,double,size_t,std::string,real_ps,real_ps>());
//   __tc_common_instantiation(tc_add_ps_to_arg);
// }

// TODO: separate abs into power function instantiation, and use type traits or overloading or specialization to
// establish if we need to provide abs too (since we don't need it for reals)?
// template <class T>
// void ps_instantiate_complex_specifics(class_<T> &complex_inst)
// {
//   typedef T complex_type;
//   typedef typename complex_type::value_type value_type;
// // Additional ctor(s) for complex series
//   complex_inst.def(init<int,int>());
//   complex_inst.def(init<const std::complex<int> &>());
//   complex_inst.def(init<const double &,const double &>());
//   complex_inst.def(init<const std::complex<double> &>());
//   complex_inst.def(init<value_type>());
//   complex_inst.def(init<value_type,value_type>());
//   complex_inst.def("real", &complex_type::real);
//   complex_inst.def("imag", &complex_type::imag);
// // TODO: this must be moved into own toolbox.
// //  complex_inst.def("abs", &complex_type::abs);
//   complex_inst.def("conj", &complex_type::conj);
// // Assignment.
// //   complex_inst.def(self=std::complex<int>());
// //   complex_inst.def(self=std::complex<double>());
// //   complex_inst.def(self=value_type());
// // Addition and subtraction.
//   complex_inst.def(self+=std::complex<int>());
//   complex_inst.def(self+=std::complex<double>());
//   complex_inst.def(self+=value_type());
//   complex_inst.def(self+std::complex<int>());
//   complex_inst.def(self+std::complex<double>());
//   complex_inst.def(self+value_type());
//   complex_inst.def(self-=std::complex<int>());
//   complex_inst.def(self-=std::complex<double>());
//   complex_inst.def(self-=value_type());
//   complex_inst.def(self-std::complex<int>());
//   complex_inst.def(self-std::complex<double>());
//   complex_inst.def(self-value_type());
// // Multiplication.
//   complex_inst.def(self*=std::complex<int>());
//   complex_inst.def(self*=std::complex<double>());
//   complex_inst.def(self*=value_type());
//   complex_inst.def(self*std::complex<int>());
//   complex_inst.def(self*std::complex<double>());
//   complex_inst.def(self*value_type());
// // Division.
//   complex_inst.def(self/=std::complex<int>());
//   complex_inst.def(self/=std::complex<double>());
//   complex_inst.def(self/std::complex<int>());
//   complex_inst.def(self/std::complex<double>());
// }

// template <class T>
// void instantiate_tass17()
// {
//   class_<tass17<T> >("tass17","TASS theory, version 1.7.",no_init)
//     .def("load", &tass17<T>::load,"Load series into memory.")
//     .staticmethod("load")
//     .def("status", &tass17<T>::status,"Print to screen TASS status.")
//     .staticmethod("status")
//     .def("add_delta_lambdas", &tass17<T>::add_delta_lambdas,"Correct series with long term perturbations.")
//     .staticmethod("add_delta_lambdas")
//     .def("m0", &tass17<T>::m0,return_value_policy<copy_const_reference>(),"Get Saturn's mass in Sun mass units.")
//     .staticmethod("m0")
//     .def("m6", &tass17<T>::m6,return_value_policy<copy_const_reference>(),"Get inverse of Titan's mass in Saturn mass units.")
//     .staticmethod("m6")
//     .def("lambda4", &tass17<T>::lambda4,return_value_policy<copy_const_reference>(),"Get Dione's lambda.")
//     .staticmethod("lambda4")
//     .def("lambda6", &tass17<T>::lambda6,return_value_policy<copy_const_reference>(),"Get Titan's lambda.")
//     .staticmethod("lambda6")
//     .def("r6", &tass17<T>::r6,"Calculate Titan's radius.")
//     .staticmethod("r6")
//     .def("p6", &tass17<T>::p6,return_value_policy<copy_const_reference>(),"Get Titan's p.")
//     .staticmethod("p6")
//     .def("z6", &tass17<T>::z6,return_value_policy<copy_const_reference>(),"Get Titan's z.")
//     .staticmethod("z6")
//     .def("zeta6", &tass17<T>::zeta6,return_value_policy<copy_const_reference>(),"Get Titan's (greek) zeta.")
//     .staticmethod("zeta6")
//     .def("e", &tass17<T>::e,"Calculate eccentricity from elliptic orbital element z.")
//     .staticmethod("e")
//     .def("a", &tass17<T>::a,"Calculate semi-major axis a from elliptic orbital element p.")
//     .staticmethod("a")
//     .def("eiM", &tass17<T>::eiM,"Calculate complex exponential of mean mean motion M.")
//     .staticmethod("eiM")
//     .def("vienne_r", &tass17<T>::vienne_r,"Calculate radius using Vienne's FORTRAN routine.")
//     .staticmethod("vienne_r");
//   class_<tc_vienne_r6<T> > tc_vienne_r6_inst("tc_vienne_r6",
//     init<typename tc_vienne_r6<T>::b_type,double,double,size_t>());
//   __tc_common_instantiation(tc_vienne_r6);
// #ifdef _PIRANHA_TBB
// // We need also the non-parallel evaluator when parallel is enabled, because otherwise we won't be able
// // to plot from Python.
//   range_evaluator_instantiation<tc_vienne_r6<T>,false>("__range_evaluator_tc_vienne_r6_non_parallel",
//     "Non-parallel range evaluator for tc_vienne_r6.");
// #endif
// }

#endif
