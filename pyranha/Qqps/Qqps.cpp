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
#include <boost/python/docstring_options.hpp>
#include <string>

#include "../../src/manipulators/qqps.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace piranha::manipulators;
using namespace pyranha;

BOOST_PYTHON_MODULE(_Qqps)
{
    docstring_options docOptions(true, false, false);
    translate_exceptions();

    class_<qqps> inst(
        series_basic_instantiation<qqps>(std::string("qqps"),
        std::string("Poisson series with arbitrary-precision rational coefficients and arbitrary-precision rational exponents.")));
    common_poisson_series_instantiation(inst, "qqps");
    celmec_instantiation(inst);
    series_trigonometric_instantiation(inst);
    series_sub_instantiation<qqps, qqps>(inst);
    series_ei_sub_instantiation<qqps, qqpsc>(inst);
    class_<qqpsc> instc(
        series_basic_instantiation<qqpsc>(std::string("qqpsc"),
        std::string("Poisson series with complex arbitrary-precision "
        "rational coefficients and arbitrary-precision rational exponents.")));
    common_poisson_series_instantiation(instc, "qqpsc");
    series_complex_instantiation(instc, inst);
    series_sub_instantiation<qqpsc, qqps>(instc);
    series_ei_sub_instantiation<qqpsc, qqpsc>(instc);
}
