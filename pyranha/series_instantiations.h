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

#ifndef PYRANHA_SERIES_INSTANTIATIONS_H
#define PYRANHA_SERIES_INSTANTIATIONS_H

#include <boost/python/class.hpp>
#include <boost/python/operators.hpp>
#include <string>

#include "../src/core/integer_typedefs.h"
#include "../src/core/psym.h"

namespace pyranha
{
  /// Basic series instantiation.
  template <class T>
    boost::python::class_<T> series_basic_instantiation(const std::string &name, const std::string &description)
  {
    boost::python::class_<T> inst(name.c_str(),description.c_str());
    inst.def(boost::python::init<const T &>());
    inst.def(boost::python::init<const std::string &>());
    inst.def(boost::python::init<const piranha::max_fast_int &>());
    inst.def(boost::python::init<const double &>());
    inst.def("__copy__",&T::copy);
    inst.def("__repr__",&T::print_to_string);
    inst.def("__len__",&T::length);
    inst.def("save_to",&T::save_to,"Save series to file.");
    inst.def("eval",&T::eval);
    inst.def("swap",&T::swap);
    // NOTICE: the order seems important here, if we place *=int before *=double we
    // will get just *=double in Python. Go figure...
    // Addition and subtraction.
    inst.def(boost::python::self+=piranha::max_fast_int());
    inst.def(boost::python::self+=double());
    inst.def(boost::python::self+=boost::python::self);
    inst.def(boost::python::self+piranha::max_fast_int());
    inst.def(piranha::max_fast_int()+boost::python::self);
    inst.def(boost::python::self+double());
    inst.def(double()+boost::python::self);
    inst.def(boost::python::self+boost::python::self);
    inst.def(boost::python::self-=piranha::max_fast_int());
    inst.def(boost::python::self-=double());
    inst.def(boost::python::self-=boost::python::self);
    inst.def(boost::python::self-piranha::max_fast_int());
    inst.def(piranha::max_fast_int()-boost::python::self);
    inst.def(boost::python::self-double());
    inst.def(double()-boost::python::self);
    inst.def(boost::python::self-boost::python::self);
    // Multiplication.
    inst.def(boost::python::self*=piranha::max_fast_int());
    inst.def(boost::python::self*=double());
    inst.def(boost::python::self*=boost::python::self);
    inst.def(boost::python::self*piranha::max_fast_int());
    inst.def(piranha::max_fast_int()*boost::python::self);
    inst.def(boost::python::self*double());
    inst.def(double()*boost::python::self);
    inst.def(boost::python::self*boost::python::self);
    // Division.
    inst.def(boost::python::self/=piranha::max_fast_int());
    inst.def(boost::python::self/=double());
    inst.def(boost::python::self/piranha::max_fast_int());
    inst.def(boost::python::self/double());
    // Exponentiation.
    inst.def("__pow__",&T::pow);
    return inst;
  }

  template <class T>
    void series_trigonometric_instantiation(boost::python::class_<T> &inst)
  {
    inst.def("cos",&T::cos);
    inst.def("sin",&T::sin);
  }

  template <class T>
    void series_differential_instantiation(boost::python::class_<T> &inst)
  {
    typedef T (T::*partial_name)(const std::string &) const;
    typedef T (T::*partial_psym)(const piranha::psym &) const;
    inst.def("partial",partial_name(&T::partial));
    inst.def("partial",partial_psym(&T::partial));
  }

  template <class T>
    void series_psym_instantiation(boost::python::class_<T> &inst)
  {
    inst.def(boost::python::init<const piranha::psym &>());
  }

  template <class T>
    void series_special_functions_instantiation(boost::python::class_<T> &inst)
  {
    inst.def("besselJ",&T::besselJ,"Bessel function of the first kind of integer order.");
    inst.def("dbesselJ",&T::dbesselJ,"Partial derivative of Bessel function of the first kind of integer order.");
  }

  template <class T>
    void celmec_instantiation(boost::python::class_<T> &inst)
  {
    inst.def("r_a",&T::r_a,"Elliptic expansion of r / a.").staticmethod("r_a");
    inst.def("sin_f",&T::sin_f,"Ellipic expansion of sin(f).").staticmethod("sin_f");
  }

  template <class T>
    void power_series_instantiation(boost::python::class_<T> &inst)
  {
    inst.def("degree",&T::degree,"Get the degree of the power series.");
    inst.def("min_degree",&T::min_degree,"Get the minimum degree of the power series.");
  }

  template <class T>
    void common_polynomial_instantiation(boost::python::class_<T> &inst)
  {
    series_psym_instantiation(inst);
    series_differential_instantiation(inst);
    series_special_functions_instantiation(inst);
    power_series_instantiation(inst);
  }

  template <class T>
    void common_poisson_series_instantiation(boost::python::class_<T> &inst)
  {
    series_trigonometric_instantiation(inst);
    series_psym_instantiation(inst);
    series_differential_instantiation(inst);
    series_special_functions_instantiation(inst);
    power_series_instantiation(inst);
    celmec_instantiation(inst);
  }

  template <class T>
    void fourier_specifics(boost::python::class_<T> &inst)
  {}
}

#endif
