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


#include "../../src/manipulators/zpoly.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

#include "pybind11/pybind11.h"


PYBIND11_MODULE(_Zpoly, m)
{
    //docstring_options docOptions(true, false, false);
    //translate_exceptions();  //still todo

    pybind11::class_<piranha::manipulators::zpoly> inst(pyranha::series_basic_instantiation<piranha::manipulators::zpoly>(m, 
                "zpoly", "Multivariate polynomial with arbitrary-size integer coefficients."));
    pyranha::common_polynomial_instantiation(inst);
    pyranha::series_sub_instantiation<piranha::manipulators::zpoly, piranha::manipulators::zpoly>(inst);

    // complex zpoly
    pybind11::class_<piranha::manipulators::zpolyc> instc(pyranha::series_basic_instantiation<piranha::manipulators::zpolyc>(m,
                "zpolyc", "Multivariate polynomial with complex arbitrary-size integer coefficients."));
    pyranha::common_polynomial_instantiation(instc);
    pyranha::series_complex_instantiation(instc, inst);
    pyranha::series_sub_instantiation<piranha::manipulators::zpolyc, piranha::manipulators::zpoly>(instc);
    pyranha::series_sub_instantiation<piranha::manipulators::zpolyc, piranha::manipulators::zpolyc>(instc);
}
