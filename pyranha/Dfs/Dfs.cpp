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
#include "../series_instantiations.h"
#include "../exceptions.h"

#include "pybind11/pybind11.h"


PYBIND11_MODULE(_Dfs, m)
{
    //docstring_options docOptions(true, false, false);
    pyranha::translate_exceptions();

    pybind11::class_<piranha::manipulators::dfs> inst(pyranha::series_basic_instantiation<piranha::manipulators::dfs>(m, "dfs", "Fourier series with double precision coefficients."));
    pyranha::common_fourier_series_instantiation(inst);
    pyranha::series_trigonometric_instantiation(inst);
    pyranha::series_sub_instantiation<piranha::manipulators::dfs, piranha::manipulators::dfs>(inst);

    pybind11::class_<piranha::manipulators::dfsc> instc(pyranha::series_basic_instantiation<piranha::manipulators::dfsc>(m, "dfsc", "Fourier series with complex double precision coefficients."));
    pyranha::series_complex_instantiation(instc, inst);
    pyranha::common_fourier_series_instantiation(instc);
    pyranha::series_sub_instantiation<piranha::manipulators::dfsc, piranha::manipulators::dfs>(instc);
}
