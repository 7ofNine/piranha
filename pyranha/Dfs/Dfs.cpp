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

#include "../../src/manipulators/dfs.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace piranha::manipulators;
using namespace pyranha;

BOOST_PYTHON_MODULE(_Dfs)
{
	translate_exceptions();

	class_<dfs> inst = series_basic_instantiation<dfs>(std::string("dfs"),
			std::string("Fourier series with double precision coefficients."));
	common_fourier_series_instantiation(inst);
	series_trigonometric_instantiation(inst);
	series_sub_instantiation<dfs, dfs>(inst);
	class_<dfsc> instc = series_basic_instantiation<dfsc>(std::string("dfsc"),
			std::string("Fourier series with complex double precision coefficients."));
	series_complex_instantiation(instc, inst);
	common_fourier_series_instantiation(instc);
	series_sub_instantiation<dfsc, dfs>(instc);
}
