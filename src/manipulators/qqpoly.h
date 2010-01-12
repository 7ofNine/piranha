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

#ifndef PIRANHA_QQPOLY_H
#define PIRANHA_QQPOLY_H

#include <complex>

#include "../core/base_classes/series_multiplier.h"
#include "../core/numerical_coefficients/mpq_cf.h"
#include "../core/polynomial_common/q_expo_array.h"
#include "../core/polynomial/polynomial.h"
#include "../core/truncators/power_series.h"

namespace piranha
{
namespace manipulators
{
	/// Manipulator of multivariate polynomials with arbitrary-size rational coefficients.
	typedef polynomial
	<
		mpq_cf,
		q_expo_array<0>::type,
		series_multiplier,
		truncators::power_series
	> qqpoly;

	typedef std::complex<qqpoly> qqpolyc;
}
}

#endif