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

#include <boost/math/special_functions/gamma.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <complex>

#include "../../src/core/math.h"
#include "../exceptions.h"

using namespace pyranha;
using namespace piranha;
using namespace boost::python;

BOOST_PYTHON_MODULE(_Math)
{
	translate_exceptions();
	def("besselJ", &besselJ, "Bessel function of the first kind of integer order.");
	def("besselI", &besselI, "Modified Bessel function of the first kind of integer order.");
	def("Pnm", &Pnm, "Associated Legendre function.");
	def("Pn", &Pn, "Legendre polynomial.");
	def("Ynm", &Ynm, "Spherical harmonic (not normalised).");
	def("factorial", &piranha::factorial, "Factorial.");
	def("double_factorial", &piranha::double_factorial, "Double factorial of non-negative integer argument.");
	typedef double (*double_gamma)(double);
	def("gamma", double_gamma(&boost::math::tgamma<double>), "Gamma function.");
	def("cs_phase", &cs_phase, "Condon-Shortley phase (-1**m).");
}
