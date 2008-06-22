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

#include "../../src/manipulators/dps.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace piranha::manipulators;
using namespace pyranha;

BOOST_PYTHON_MODULE(_Dps)
{
	translate_exceptions();

	std::pair<class_<dps>,class_<dps::term_type> > inst = series_basic_instantiation<dps>(std::string("dps"),
									 std::string("Poisson series with double precision coefficients."));
	common_poisson_series_instantiation(inst.first);
	celmec_instantiation(inst.first);
	series_trigonometric_instantiation(inst.first);
	std::pair<class_<dpsc>,class_<dpsc::term_type> > instc = series_basic_instantiation<dpsc>(std::string("dpsc"),
									   std::string("Poisson series with complex double precision coefficients."));
	common_poisson_series_instantiation(instc.first);
	series_complex_instantiation(instc.first, inst.first);
}
