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

#ifndef PIRANHA_QPOLY_H
#define PIRANHA_QPOLY_H

#include <complex>

#include "../core/numerical_coefficients/mpq_cf.h"
#include "../core/polynomial_common/expo_vector.h"
#include "../core/polynomial/polynomial.h"
#include "../core/polynomial_common/polynomial_multiplier.h"
#include "../core/truncators/power_series.h"

namespace piranha
{
namespace manipulators
{
	/// Manipulator of multivariate polynomials with arbitrary-size rational coefficients.
	typedef Polynomial
	<
		mpq_cf,
		ExpoVector<int16_t, 0>,
		polynomial_multiplier,
		truncators::PowerSeries
	> qpoly;

	typedef std::complex<qpoly> qpolyc;
}
}

#endif
