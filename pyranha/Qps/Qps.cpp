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


PYBIND11_MODULE(_Qps, m)
{
    //docstring_options docOptions(true, false, false);
    pyranha::translate_exceptions();

    pybind11::class_<piranha::manipulators::qps> inst(pyranha::series_basic_instantiation<piranha::manipulators::qps>(m, "qps",
            "Poisson series with arbitrary-precision rational coefficients."));
    pyranha::common_poisson_series_instantiation(inst, "qps");
    pyranha::celmec_instantiation(inst);
    pyranha::series_trigonometric_instantiation(inst);
    pyranha::series_sub_instantiation<piranha::manipulators::qps, piranha::manipulators::qps>(inst);
    pyranha::series_ei_sub_instantiation<piranha::manipulators::qps, piranha::manipulators::qpsc>(inst);

    pybind11::class_<piranha::manipulators::qpsc> instc(pyranha::series_basic_instantiation<piranha::manipulators::qpsc>(m, "qpsc", "Poisson series with complex arbitrary-precision rational coefficients."));
    pyranha::common_poisson_series_instantiation(instc, "qpsc");
    pyranha::series_complex_instantiation(instc, inst);
    pyranha::series_sub_instantiation<piranha::manipulators::qpsc, piranha::manipulators::qps>(instc);
    pyranha::series_ei_sub_instantiation<piranha::manipulators::qpsc, piranha::manipulators::qpsc>(instc);
}
