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

#ifndef PIRANHA_TRIG_ARRAY_COMMONS_H
#define PIRANHA_TRIG_ARRAY_COMMONS_H

#include <boost/algorithm/string/split.hpp>
#include <cmath> // For std::abs.
#include <complex>
#include <string>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/series_builders.h"
#include "../config.h"
#include "../psym.h"
#include "../settings.h"
#include "../utils.h" // For is_integer().
#include "trig_evaluator.h"

#define derived_const_cast (static_cast<Derived const *>(this))
#define derived_cast (static_cast<Derived *>(this))

namespace piranha
{
	/// Common class for dense trigonometric array.
	/**
	 * Intended to extend piranha::int_array for the manipulation of trigonometric
	 * parts in Poisson series.
	 */
	template <class Derived>
	class trig_array_commons
	{
		public:


	};
}

#undef derived_const_cast
#undef derived_cast

#endif
