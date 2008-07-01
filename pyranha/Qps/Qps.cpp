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
#include <utility>

#include "../../src/manipulators/qps.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace piranha::manipulators;
using namespace pyranha;

BOOST_PYTHON_MODULE(_Qps)
{
	translate_exceptions();

	std::pair<class_<qps>,class_<qps::term_type> > inst(
			series_basic_instantiation<qps>(std::string("qps"),
			std::string("Poisson series with arbitrary-precision rational coefficients.")));
	common_poisson_series_instantiation(inst.first,"qps");
	celmec_instantiation(inst.first);
	series_trigonometric_instantiation(inst.first);

	std::pair<class_<qpsc>,class_<qpsc::term_type> > instc(
				series_basic_instantiation<qpsc>(std::string("qpsc"),
				std::string("Poisson series with complex arbitrary-precision rational coefficients.")));
	common_poisson_series_instantiation(instc.first,"qpsc");
	series_complex_instantiation(instc.first, inst.first);
}
