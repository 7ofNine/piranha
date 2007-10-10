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

#ifndef PYRANHA_H
#define PYRANHA_H

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/self.hpp>

#include "../src/piranha.h"

using namespace boost::python;
using namespace piranha;

/// Instantiation of spectral comparison.
template <class T>
class_<sc<T> > sc_instatiation(const std::string &name)
{
  class_<sc<T> > retval("sc",init<T,T,bool>());
  retval.def(init<T,T>());
  retval.def("size",&sc<T>::size);
  retval.def("is_relative",&sc<T>::is_relative);
  retval.def("diffs",&sc<T>::diffs,return_value_policy<copy_const_reference>());
  retval.def("gnuplot_save", &sc<T>::gnuplot_save);
  return retval;
}


/// Instantiation of common functions for time comparison.
// Here T is the tc of operation<series_type>.
template <class T>
void tc_common_instantiation(class_<T> &time_c)
{
  time_c.def("sigma",&T::sigma,return_value_policy<copy_const_reference>());
  time_c.def("max_error",&T::max_error,return_value_policy<copy_const_reference>());
  time_c.def("time",&T::time,return_value_policy<copy_const_reference>());
  time_c.def("hs",&T::hs,return_value_policy<copy_const_reference>());
  time_c.def("hs_computed",&T::hs_computed,return_value_policy<copy_const_reference>());
  time_c.def("error",&T::error);
  time_c.def("gnuplot_save",&T::gnuplot_save);
  time_c.def("size",&T::size);
}


/// Template for the instantiation of a ps. It will expose methods common to real and complex
// ps, the specializations take place below.
template <class T>
class_<T> ps_basic_instantiation(const std::string &name, const std::string &description)
{
// FIXME: we need to define all the operations on ints, double, rationals (?) etc etc.
// FIXME: ... as well as the ctors.
// This is a trick to help resolve overloaded methods inside classes.
  typedef void (T::*crop_it)(const typename T::it_s_index &);
  typedef void (T::*crop_real)(const double &);
  typedef void (T::*put_noargs)() const;
  typedef void (T::*put_n)(int) const;
  typedef void (T::*put_terms_noargs)() const;
  typedef void (T::*put_terms_n)(int) const;
  typedef void (T::*put_phases_freqs_noargs)() const;
  typedef void (T::*put_phases_freqs_n)(int) const;
  typedef typename T::eval_type (T::*mean_def) (const double &, const double &) const;
  typedef typename T::eval_type (T::*mean_n)(const double &, const double &,
    const size_t &) const;
  class_<T> inst(name.c_str(),description.c_str());
  inst.def(init<const T &>());
  inst.def(init<const std::string &>());
  inst.def(init<const double &>());
  inst.def(init<int>());
  inst.def("__copy__", &T::copy);
  inst.def("__iter__", iterator<T,return_internal_reference<> >());
  inst.def("__len__", &T::length);
  inst.def("begin", &T::begin);
  inst.def("end", &T::end);
  inst.def("address", &T::address);
  inst.def("save_to", &T::save_to, "Save series to file.");
  inst.def("put", put_noargs(&T::put));
  inst.def("put", put_n(&T::put));
  inst.def("put_terms", put_terms_noargs(&T::put_terms));
  inst.def("put_terms", put_terms_n(&T::put_terms));
  inst.def("put_phases_freqs", put_phases_freqs_noargs(&T::put_phases_freqs));
  inst.def("put_phases_freqs", put_phases_freqs_n(&T::put_phases_freqs));
  inst.def("trig_density", &T::trig_density);
  inst.def("length", &T::length);
  inst.def("trig_width", &T::trig_width);
  inst.def("g_norm", &T::g_norm);
  inst.def("footprint", &T::footprint);
  inst.def("checkup", &T::checkup);
  inst.def("discontinuity", &T::discontinuity);
  inst.def("crop", crop_real(&T::crop));
  inst.def("crop", crop_it(&T::crop));
  inst.def("spectral_cutoff", &T::spectral_cutoff);
  inst.def("cumulative_crop", &T::cumulative_crop);
  inst.def("insert_phases", &T::insert_phases);
  inst.def("add_cf_arg", &T::add_cf_arg);
  inst.def("add_trig_arg", &T::add_trig_arg);
  inst.def("t_eval", &T::t_eval);
  inst.def("t_eval_brute", &T::t_eval_brute);
  inst.def("mean", mean_def(&T::mean));
  inst.def("mean", mean_n(&T::mean));
  inst.def("swap", &T::swap);
  inst.def(self+=self);
  inst.def(self+self);
  inst.def(self+=double());
  inst.def(self+double());
  inst.def(self-=self);
  inst.def(self-self);
  inst.def(self-=double());
  inst.def(self-double());
  inst.def(self*=self);
  inst.def(self*self);
// NOTICE: the order seems important here, if we place *=int before *=double we
// will get just *=double in Python. Go figure...
  inst.def(self*=double());
  inst.def(self*double());
  inst.def(self*=int());
  inst.def(self*int());
  inst.def(self/=double());
  inst.def(self/double());
  inst.def(self/=int());
  inst.def(self/int());

// Instantiate spectral comparison.
  sc_instatiation<T>(name);
// Instantiate common time comparisons.
  class_<tc_equal<T> > tc_equal_inst("tc_equal",
    init<typename tc_equal<T>::b_type,double,double,size_t,T>());
  tc_common_instantiation(tc_equal_inst);
  class_<tc_mult<T> > tc_mult_inst("tc_mult",
    init<typename tc_mult<T>::b_type,double,double,size_t,T,T>());
  tc_common_instantiation(tc_mult_inst);
  class_<tc_insert_phases<T> > tc_insert_phases_inst("tc_insert_phases",
    init<typename tc_insert_phases<T>::b_type,double,double,size_t,phase_list,T>());
  tc_common_instantiation(tc_insert_phases_inst);

  return inst;
}

template <class T>
void ps_instantiate_differential_specifics(class_<T> &inst)
{
  inst.def("partial", &T::partial);
}

template <class T>
void ps_instantiate_real_specifics(class_<T> &real)
{
  typedef T real_ps;
  typedef void (real_ps::*real_add_ps_to_arg_index)(trig_size_t, const real_ps &);
  typedef void (real_ps::*real_add_ps_to_arg_string)(const std::string &, const real_ps &);
  real.def(init<const psymbol &, psymbol::type>());
  real.def("complexp", &real_ps::complexp);
  real.def("cosine", &real_ps::cosine);
  real.def("sine", &real_ps::sine);
  real.def("pow", &real_ps::pow);
  real.def("add_ps_to_arg", real_add_ps_to_arg_index(&real_ps::add_ps_to_arg));
  real.def("add_ps_to_arg", real_add_ps_to_arg_string(&real_ps::add_ps_to_arg));
// External functions.
//   def("kep_cosE",&astro::kep_cosE<real_ps>,"Solve Kepler's equation for cosE.");
// // NOTE: which functions does it make sense to keep here?
//   def("Pnm",&math::Pnm<real_ps>,"Legendre function of the first kind - Pnm(cos(theta)).");
//   def("Ynm",&math::Ynm<real_ps>,"Non-normalized spherical harmonic.");
//   def("wig_rot",&math::wig_rot<real_ps>,"Wigner rotation theorem for spherical harmonics.");
  class_<tc_complexp<real_ps> > tc_complexp_inst("tc_complexp",
    init<typename tc_complexp<real_ps>::b_type,double,double,size_t,real_ps>());
  tc_common_instantiation(tc_complexp_inst);
  class_<tc_cosine<real_ps> > tc_cosine_inst("tc_cosine",
    init<typename tc_cosine<real_ps>::b_type,double,double,size_t,real_ps>());
  tc_common_instantiation(tc_cosine_inst);
  class_<tc_sine<real_ps> > tc_sine_inst("tc_sine",
    init<typename tc_sine<real_ps>::b_type,double,double,size_t,real_ps>());
  tc_common_instantiation(tc_sine_inst);
  class_<tc_Pnm<real_ps> > tc_Pnm_inst("tc_Pnm",
    init<typename tc_Pnm<real_ps>::b_type,double,double,size_t,int,int,real_ps>());
  tc_common_instantiation(tc_Pnm_inst);
  class_<tc_Ynm<real_ps> > tc_Ynm_inst("tc_Ynm",
    init<typename tc_Ynm<real_ps>::b_type,double,double,size_t,int,int,real_ps,real_ps>());
  tc_common_instantiation(tc_Ynm_inst);
  class_<tc_wig_rot<real_ps> > tc_wig_rot_inst("tc_wig_rot",
    init<typename tc_wig_rot<real_ps>::b_type,double,double,size_t,int,int,real_ps,real_ps,
    real_ps,real_ps,real_ps>());
  tc_common_instantiation(tc_wig_rot_inst);
  class_<tc_pow<real_ps> > tc_pow_inst("tc_pow",
    init<typename tc_pow<real_ps>::b_type,double,double,size_t,double,real_ps>());
  tc_common_instantiation(tc_pow_inst);
  class_<tc_add_ps_to_arg<real_ps> > tc_add_ps_to_arg_inst("tc_add_ps_to_arg",
    init<typename tc_add_ps_to_arg<real_ps>::b_type,double,double,size_t,std::string,real_ps,real_ps>());
  tc_common_instantiation(tc_add_ps_to_arg_inst);
}


template <class T>
void ps_instantiate_complex_specifics(class_<T> &complex)
{
  typedef T complex_ps;
  typedef typename complex_ps::real_ps real_ps;
  complex.def(init<const complex_double &>());
  complex.def("real", &complex_ps::real);
  complex.def("imag", &complex_ps::imag);
  complex.def("abs", &complex_ps::abs);
  complex.def("conj", &complex_ps::conj);
  complex.def("make_conj", &complex_ps::make_conj);
  complex.def(self+=real_ps());
  complex.def(self+real_ps());
  complex.def(self-=real_ps());
  complex.def(self-real_ps());
  complex.def(self*=real_ps());
  complex.def(self*real_ps());
// Additional ctor(s) for complex series
  complex.def(init<real_ps>());
  complex.def(init<real_ps,real_ps>());
// FIXME: this needs to be fixed.
  //complex.def(init<typename real_ps::cf_type,typename real_ps::cf_type>());
}

template <class T>
void instantiate_tass17()
{
  class_<tass17<T> >("tass17","TASS theory, version 1.7.",no_init)
    .def("load", &tass17<T>::load,"Load series into memory.")
    .staticmethod("load")
    .def("status", &tass17<T>::status,"Print to screen TASS status.")
    .staticmethod("status")
    .def("add_delta_lambdas", &tass17<T>::add_delta_lambdas,"Correct series with long term perturbations.")
    .staticmethod("add_delta_lambdas")
    .def("m0", &tass17<T>::m0,return_value_policy<copy_const_reference>(),"Get Saturn's mass in Sun mass units.")
    .staticmethod("m0")
    .def("m6", &tass17<T>::m6,return_value_policy<copy_const_reference>(),"Get inverse of Titan's mass in Saturn mass units.")
    .staticmethod("m6")
    .def("lambda4", &tass17<T>::lambda4,return_value_policy<copy_const_reference>(),"Get Dione's lambda.")
    .staticmethod("lambda4")
    .def("lambda6", &tass17<T>::lambda6,return_value_policy<copy_const_reference>(),"Get Titan's lambda.")
    .staticmethod("lambda6")
    .def("r6", &tass17<T>::r6,"Calculate Titan's radius.")
    .staticmethod("r6")
    .def("p6", &tass17<T>::p6,return_value_policy<copy_const_reference>(),"Get Titan's p.")
    .staticmethod("p6")
    .def("z6", &tass17<T>::z6,return_value_policy<copy_const_reference>(),"Get Titan's z.")
    .staticmethod("z6")
    .def("zeta6", &tass17<T>::zeta6,return_value_policy<copy_const_reference>(),"Get Titan's (greek) zeta.")
    .staticmethod("zeta6")
    .def("e", &tass17<T>::e,"Calculate eccentricity from elliptic orbital element z.")
    .staticmethod("e")
    .def("a", &tass17<T>::a,"Calculate semi-major axis a from elliptic orbital element p.")
    .staticmethod("a")
    .def("eiM", &tass17<T>::eiM,"Calculate complex exponential of mean mean motion M.")
    .staticmethod("eiM")
    .def("vienne_r", &tass17<T>::vienne_r,"Calculate radius using Vienne's FORTRAN routine.")
    .staticmethod("vienne_r")
    ;

  class_ <tc_vienne_r6<T> >
    ((std::string("tc_vienne_r6")).c_str(),
    init<typename tc_vienne_r6<T>::b_type,double,double,size_t>())
    .def("sigma",&tc_vienne_r6<T>::sigma,return_value_policy<copy_const_reference>())
    .def("max_error",&tc_vienne_r6<T>::max_error,
    return_value_policy<copy_const_reference>())
    .def("time",&tc_vienne_r6<T>::time,return_value_policy<copy_const_reference>())
    .def("hs",&tc_vienne_r6<T>::hs,return_value_policy<copy_const_reference>())
    .def("hs_computed",&tc_vienne_r6<T>::hs_computed,
    return_value_policy<copy_const_reference>())
    .def("error",&tc_vienne_r6<T>::error)
    .def("gnuplot_save",&tc_vienne_r6<T>::gnuplot_save)
    .def("size",&tc_vienne_r6<T>::size)
    ;
}

#endif
