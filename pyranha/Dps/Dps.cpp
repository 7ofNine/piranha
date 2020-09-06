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


#include "../../src/manipulators/dfs.h"
#include "../../src/manipulators/dps.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

#include"pybind11/pybind11.h"


PYBIND11_MODULE(_Dps, m)
{
    //docstring_options docOptions(true, false, false);
    pyranha::translate_exceptions();

    pybind11::class_<piranha::manipulators::dps> inst(pyranha::series_basic_instantiation<piranha::manipulators::dps>(m, "dps", "Poisson series with double precision coefficients."));
    pyranha::common_poisson_series_instantiation(inst, "dps");
    pyranha::celmec_instantiation(inst);
    pyranha::series_trigonometric_instantiation(inst);
    pyranha::series_sub_instantiation<piranha::manipulators::dps, piranha::manipulators::dps>(inst);
    pyranha::series_ei_sub_instantiation<piranha::manipulators::dps, piranha::manipulators::dpsc>(inst);

    inst.def("to_dfs", &piranha::manipulators::dps::to_fs<piranha::manipulators::dfs>, "Convert to dfs.");

    pybind11::class_<piranha::manipulators::dpsc> instc(pyranha::series_basic_instantiation<piranha::manipulators::dpsc>(m, "dpsc", "Poisson series with complex double precision coefficients."));
    pyranha::common_poisson_series_instantiation(instc, "dpsc");
    pyranha::series_complex_instantiation(instc, inst);
    pyranha::series_sub_instantiation<piranha::manipulators::dpsc, piranha::manipulators::dps>(instc);
    pyranha::series_ei_sub_instantiation<piranha::manipulators::dpsc, piranha::manipulators::dpsc>(instc);
    instc.def("to_dfsc", &piranha::manipulators::dpsc::to_fs<piranha::manipulators::dfsc>, "Convert to dfsc.");
}
