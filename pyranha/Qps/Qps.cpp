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

#include "../../src/manipulators/qps.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

#include "pybind11/pybind11.h"

using namespace piranha;
using namespace piranha::manipulators;
using namespace pyranha;

PYBIND11_MODULE(_Qps, m)
{
    //docstring_options docOptions(true, false, false);
    //translate_exceptions();

    pybind11::class_<qps> inst(series_basic_instantiation<qps>(m, "qps",
            "Poisson series with arbitrary-precision rational coefficients."));
    common_poisson_series_instantiation(inst, "qps");
    celmec_instantiation(inst);
    series_trigonometric_instantiation(inst);
    series_sub_instantiation<qps, qps>(inst);
    series_ei_sub_instantiation<qps, qpsc>(inst);

    pybind11::class_<qpsc> instc(series_basic_instantiation<qpsc>(m, "qpsc", "Poisson series with complex arbitrary-precision rational coefficients."));
    common_poisson_series_instantiation(instc, "qpsc");
    series_complex_instantiation(instc, inst);
    series_sub_instantiation<qpsc, qps>(instc);
    series_ei_sub_instantiation<qpsc, qpsc>(instc);
}
