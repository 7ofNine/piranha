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

#include <boost/python/class.hpp>
#include <boost/python/operators.hpp>
#include <string>

#include "commons.h"

// TODO
// - expose complex mp types, including choose overloads for complexes.

namespace pyranha
{
	template <class T, class U>
	inline void common_mp_operators(T &inst, const U &x)
	{
		// Comparisons
		inst.def(boost::python::self == x);
		inst.def(boost::python::self != x);
		inst.def(x == boost::python::self);
		inst.def(x != boost::python::self);

		// Addition and subtraction.
		inst.def(boost::python::self += x);
		inst.def(boost::python::self + x);
		inst.def(x + boost::python::self);
		inst.def(boost::python::self -= x);
		inst.def(boost::python::self - x);
		inst.def(x - boost::python::self);

		// Multiplication.
		inst.def(boost::python::self *= x);
		inst.def(boost::python::self * x);
		inst.def(x * boost::python::self);

		// Division.
		inst.def(boost::python::self /= x);
		inst.def(boost::python::self / x);
		inst.def(x / boost::python::self);
	}

	template <class T>
	inline void common_mp_methods(boost::python::class_<T>  &inst)
	{
		// Constructors.
		inst.def(boost::python::init<const int &>());
		inst.def(boost::python::init<const double &>());
		inst.def(boost::python::init<const std::string &>());
		inst.def(boost::python::init<const piranha::mp_rational &>());
		inst.def(boost::python::init<const piranha::mp_integer &>());

		// Some special methods.
		inst.def("__abs__",  &T::abs, "Absolute value.");
		inst.def("__copy__", &py_copy<T>);
		inst.def(boost::python::self_ns::repr(boost::python::self));

		// Negation.
		inst.def(-boost::python::self);

		// Exponentiation & friends.
		typedef T (T::*pow_double)(const double &) const;
		typedef T (T::*pow_rational)(const piranha::mp_rational &) const;
		inst.def("__pow__", (pow_double)&T::pow,        "Exponentiation.");
		inst.def("__pow__", (pow_rational)&T::pow,      "Exponentiation.");
		inst.def("_latex_", &py_print_to_string_tex<T>, "Latex representation.");
		inst.def("root",    &T::root,                   "N-th root.");
		inst.def("__hash__",&T::hash,                   "Hash value.");
	}

	template <class T, class U>
	inline void real_mp_operators(T &inst, const U &x)
	{
		common_mp_operators(inst,x);
		// Comparisons
		inst.def(boost::python::self < x);
		inst.def(boost::python::self <= x);
		inst.def(boost::python::self > x);
		inst.def(boost::python::self >= x);
		inst.def(x < boost::python::self);
		inst.def(x <= boost::python::self);
		inst.def(x > boost::python::self);
		inst.def(x >= boost::python::self);
	}

	template <class T>
	inline void real_mp_methods(boost::python::class_<T> &inst)
	{
		common_mp_methods(inst);
		// Conversions.
		inst.def("__float__", &T::to_double, "Convert to floating point.");
		inst.def("__int__",   &T::to_int,    "Convert to integer.");
	}

	template <class T>
	inline boost::python::class_<T> expose_real_mp_class(const std::string &name, const std::string &doc)
	{
		boost::python::class_<T> inst(name.c_str(),doc.c_str(),boost::python::init<>());
		real_mp_methods(inst);
		// Operators against standard types.
		real_mp_operators(inst,double());
		real_mp_operators(inst,int());
		real_mp_operators(inst,piranha::mp_rational());
		real_mp_operators(inst,piranha::mp_integer());
		return inst;
	}
}

#endif
