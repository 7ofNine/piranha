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

#ifndef PIRANHA_DQPS_H
#define PIRANHA_DQPS_H

#include <boost/cstdint.hpp>
#include <complex>

#include "../core/mp.h"
#include "../core/harmonic_series/trig_vector.h"
#include "../core/numerical_coefficients/double_cf.h"
#include "../core/polynomial_common/expo_vector.h"
#include "../core/poisson_series_common/poisson_series_multiplier.h"
#include "../core/poisson_series/poisson_series.h"
#include "../core/polynomial_common/polynomial_multiplier.h"
#include "../core/truncators/power_series.h"

namespace piranha
{
namespace manipulators
{
	typedef poisson_series
	<
		double_cf,
		ExpoVector<mp_rational,0>,
		TrigVector<boost::int16_t, 1>,
		polynomial_multiplier,
		poisson_series_multiplier,
		truncators::PowerSeries,
		truncators::PowerSeries
	> dqps;

	typedef std::complex<dqps> dqpsc;
}
}

#endif
