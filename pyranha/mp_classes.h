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

namespace pyranha
{
	template <class T>
	inline boost::python::class_<T> expose_mp_class(const std::string &name, const std::string &doc)
	{
		boost::python::class_<T> inst(name.c_str(),doc.c_str(),boost::python::init<>());
		inst.def(boost::python::init<const int &>());
		inst.def(boost::python::init<const double &>());
		inst.def(boost::python::init<const std::string &>());
		// Some special methods.
		inst.def("__copy__", &py_copy<T>);
		inst.def(boost::python::self_ns::repr(boost::python::self));
		inst.def(boost::python::self_ns::float_(boost::python::self));
		inst.def(boost::python::self_ns::int_(boost::python::self));
		// Equality.
		inst.def(boost::python::self == int());
		inst.def(boost::python::self == double());
		inst.def(boost::python::self == boost::python::self);
		inst.def(boost::python::self != int());
		inst.def(boost::python::self != double());
		inst.def(boost::python::self != boost::python::self);
		// Addition and subtraction.
		inst.def(boost::python::self += int());
		inst.def(boost::python::self += double());
		inst.def(boost::python::self += boost::python::self);
		inst.def(boost::python::self + int());
		inst.def(boost::python::self + double());
		inst.def(int() + boost::python::self);
		inst.def(double() + boost::python::self);
		inst.def(boost::python::self + boost::python::self);
		inst.def(boost::python::self -= int());
		inst.def(boost::python::self -= double());
		inst.def(boost::python::self -= boost::python::self);
		inst.def(boost::python::self - int());
		inst.def(boost::python::self - double());
		inst.def(int() - boost::python::self);
		inst.def(double() - boost::python::self);
		inst.def(boost::python::self - boost::python::self);
		inst.def(-boost::python::self);
		// Multiplication.
		inst.def(boost::python::self *= int());
		inst.def(boost::python::self *= double());
		inst.def(boost::python::self *= boost::python::self);
		inst.def(boost::python::self * int());
		inst.def(boost::python::self * double());
		inst.def(int() * boost::python::self);
		inst.def(double() * boost::python::self);
		inst.def(boost::python::self * boost::python::self);
		// Division.
		inst.def(boost::python::self /= int());
		inst.def(boost::python::self /= double());
		inst.def(boost::python::self /= boost::python::self);
		inst.def(boost::python::self / int());
		inst.def(boost::python::self / double());
		inst.def(boost::python::self / boost::python::self);
		// Exponentiation & friends.
		inst.def("__pow__",&T::pow,"Exponentiation.");
		inst.def("root",&T::root,"N-th root.");
		return inst;
	}
}

#endif
