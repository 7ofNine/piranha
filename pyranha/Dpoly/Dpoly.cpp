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
#include <boost/python/module.hpp>
#include <string>

#include "../../src/manipulators/dpoly.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

using namespace pyranha;
using namespace piranha;
using namespace boost::python;

BOOST_PYTHON_MODULE(_Dpoly)
{
	translate_exceptions();

	class_<manipulators::dpoly> inst = series_basic_instantiation<manipulators::dpoly>(std::string("dpoly"),
									   std::string("Multivariate polynomial with double precision coefficients."));
	//series_trigonometric_instantiation(inst);
	common_polynomial_instantiation(inst);
	class_<manipulators::dpolyc> instc = series_basic_instantiation<manipulators::dpolyc>(std::string("dpolyc"),
									  std::string("Multivariate polynomial with complex double precision coefficients."));
	series_complex_instantiation(instc);
	//ps_instantiate_differential_specifics(inst);
	/*ps_instantiate_real_specifics(inst);
	def("pow_besselJ",math::pow_besselJ<gsp,mpz_class>,
	  "Bessel function of the first kind, power series implementation.");*/
}
