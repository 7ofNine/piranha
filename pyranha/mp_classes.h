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

#ifndef PYRANHA_MP_CLASSES_H
#define PYRANHA_MP_CLASSES_H

//#include <boost/python/class.hpp>
//#include <boost/python/operators.hpp>
#include <string>

#include "commons.h"


#include "pybind11/pybind11.h"
#include "pybind11/operators.h"

// TODO
// - expose complex mp types, including choose overloads for complexes.

namespace pyranha
{
	template <class T, class U>
	inline void common_mp_operators(T &inst, const U &x)
	{
		// Comparisons
		inst.def(pybind11::self == x);
		inst.def(pybind11::self != x);
		inst.def(x == pybind11::self);   
		inst.def(x != pybind11::self);

		//// Addition and subtraction.
		inst.def(pybind11::self += x);
		inst.def(pybind11::self + x);
		inst.def(x + pybind11::self);
		inst.def(pybind11::self -= x);
		inst.def(pybind11::self - x);
		inst.def(x - pybind11::self);

		//// Multiplication.
		inst.def(pybind11::self *= x);
		inst.def(pybind11::self * x);
		inst.def(x * pybind11::self);

		//// Division.
		inst.def(pybind11::self /= x);
		inst.def(pybind11::self / x);
		inst.def(x / pybind11::self);
	}


	template <class T>
	inline void common_mp_methods(pybind11::class_<T>  &inst)
	{
		// Constructors.
		inst.def(pybind11::init<const int &>());
		inst.def(pybind11::init<const double &>());
		inst.def(pybind11::init<const std::string &>());
		//inst.def(pybind11::init<const piranha::mp_rational &>());/// how do I actually instantiat this way?? 
		//inst.def(pybind11::init<const piranha::mp_integer &>()); /// how do I actually instantiat this way??

		inst.def("__repr__", &py_print_to_string<T>, "Display function");  // make it displayable

		// Some special methods.
		inst.def("__abs__",  &T::abs, "Absolute value."); 
		inst.def("__copy__", &py_copy<T>);
		//inst.def(pybind11::self_ns::repr(pybind11::self)); //TODO: what might correspond to self_ns in pybind what does that do anyways

		// Negation.
		inst.def(- pybind11::self);

		// Exponentiation & friends.
		typedef T (T::*pow_double)(const double &) const;
		//typedef T (T::*pow_rational)(const piranha::mp_rational &) const;
		inst.def("__pow__", (pow_double)&T::pow,        "Exponentiation.");
		//inst.def("__pow__", (pow_rational)&T::pow,      "Exponentiation.");
		inst.def("_latex_", &py_print_to_string_tex<T>, "Latex representation.");
		inst.def("root",    &T::root,                   "N-th root.");
		inst.def("__hash__",&T::hash,                   "Hash value.");
	}


	template <class T, class U>
	inline void real_mp_operators(T &inst, const U &x)
	{
		common_mp_operators(inst, x);
		// Comparisons
		//inst.def(pybind11::self < x);
		//inst.def(pybind11::self <= x);
		//inst.def(pybind11::self > x);
		//inst.def(pybind11::self >= x);
		//inst.def(x < pybind11::self);
		//inst.def(x <= pybind11::self);
		//inst.def(x > pybind11::self);
		//inst.def(x >= pybind11::self);
	}


	template <class T>
	inline void real_mp_methods(pybind11::class_<T> &inst)
	{
		common_mp_methods(inst);
		// Conversions.
		inst.def("__float__", &T::to_double,  "Convert to floating point.");
		inst.def("__int__",   &T::to_int,     "Convert to integer.");
	}


	template <class T>
	inline pybind11::class_<T> expose_real_mp_class(pybind11::module m, const std::string &name, const std::string &doc)
	{
		pybind11::class_<T> inst(m, name.c_str(), doc.c_str());
	
		inst.def(pybind11::init<>());

		real_mp_methods(inst);

		// Operators against standard types.
		real_mp_operators(inst, double());
		real_mp_operators(inst, int());
		//real_mp_operators(inst, piranha::mp_rational());   // How do we create this case???
		//real_mp_operators(inst, piranha::mp_integer());    // how dowe create this case??  some method that returns such a result??
		return inst;
	}
}

#endif
