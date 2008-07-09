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

#include "../../src/manipulators/zpoly.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace piranha::manipulators;
using namespace pyranha;

BOOST_PYTHON_MODULE(_Zpoly)
{
	translate_exceptions();

	std::pair<class_<zpoly>,class_<zpoly::term_type> > inst = series_basic_instantiation<zpoly>(std::string("zpoly"),
									   std::string("Multivariate polynomial with arbitrary-size integer coefficients."));
	common_polynomial_instantiation(inst.first);
	series_sub_instantiation<zpoly,zpoly>(inst.first);
	std::pair<class_<zpolyc>,class_<zpolyc::term_type> > instc = series_basic_instantiation<zpolyc>(std::string("zpolyc"),
									  std::string("Multivariate polynomial with complex arbitrary-size integer coefficients."));
	common_polynomial_instantiation(instc.first);
	series_complex_instantiation(instc.first, inst.first);
	series_sub_instantiation<zpolyc,zpoly>(instc.first);
}
