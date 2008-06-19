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

#ifndef PIRANHA_DFS_H
#define PIRANHA_DFS_H

#include "../core/base_classes/norm_truncator.h"
#include "../core/fourier_series/fourier_series.h"
#include "../core/numerical_coefficients/double_cf.h"
#include "../core/poisson_series_common/poisson_series_multiplier.h"
#include "../core/poisson_series_common/trig_array.h"

namespace piranha
{
	namespace manipulators {
		/// Fourier series manipulator.
		typedef fourier_series
		<
		double_cf, trig_array<16, 0>,
		poisson_series_multiplier,
		norm_truncator
		> dfs;

		typedef std::complex<dfs> dfsc;
	}
}

#endif
