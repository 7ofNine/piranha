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


#include "../../src/manipulators/dqps.h"
#include "../series_instantiations.h"
#include "../exceptions.h"

#include "pybind11/pybind11.h"


PYBIND11_MODULE(_Dqps, m)
{
    //docstring_options docOptions(true, false, false);
	pyranha::translate_exceptions();

	pybind11::class_<piranha::manipulators::dqps> inst(pyranha::series_basic_instantiation<piranha::manipulators::dqps>(m, "dqps",
					"Poisson series with double precision rational coefficients and arbitrary-precision rational exponents."));
	pyranha::common_poisson_series_instantiation(inst, "dqps");
	pyranha::celmec_instantiation(inst);
	pyranha::series_trigonometric_instantiation(inst);
	pyranha::series_sub_instantiation<piranha::manipulators::dqps, piranha::manipulators::dqps>(inst);
	pyranha::series_ei_sub_instantiation<piranha::manipulators::dqps, piranha::manipulators::dqpsc>(inst);

	pybind11::class_<piranha::manipulators::dqpsc> instc(pyranha::series_basic_instantiation<piranha::manipulators::dqpsc>(m, "dqpsc", "Poisson series with complex double precision "
		"rational coefficients and arbitrary-precision rational exponents."));
	pyranha::common_poisson_series_instantiation(instc, "dqpsc");
	pyranha::series_complex_instantiation(instc, inst);
	pyranha::series_sub_instantiation<piranha::manipulators::dqpsc, piranha::manipulators::dqps>(instc);
	pyranha::series_ei_sub_instantiation<piranha::manipulators::dqpsc, piranha::manipulators::dqpsc>(instc);
}
