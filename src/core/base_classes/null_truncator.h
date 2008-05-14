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

#ifndef PIRANHA_NULL_TRUNCATOR_H
#define PIRANHA_NULL_TRUNCATOR_H

#include "../exceptions.h"

namespace piranha
{
	/// Truncator which does not truncate.
	template <class Multiplier>
	struct null_truncator {
		null_truncator(const Multiplier &) {}
		template <class Result>
		bool accept(const Result &) const {
			return true;
		}
		template <class Cf1, class Cf2, class Key>
		bool skip(const Cf1 &, const Key &, const Cf2 &, const Key &) const {
			return false;
		}
		// Limit of a power series development of a power series.
		template <class PowerSeries, class ArgsTuple>
		static size_t power_series_limit(const PowerSeries &, const ArgsTuple &,
										 const int &start = 0, const int &step_size = 1) {
			(void)start;
			(void)step_size;
			throw unsuitable("Null truncator cannot provide limit for power series.");
		}
	};
}

#endif