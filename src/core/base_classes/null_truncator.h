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
	class null_truncator
	{
		public:
			template <class Multiplier>
			class get_type
			{
				public:
					typedef get_type type;
					get_type(const Multiplier &) {}
					template <class Result>
					bool accept(const Result &) const {
						return true;
					}
					template <class Term1, class Term2>
					bool skip(const Term1 &, const Term2 &) const {
						return false;
					}
					// Limit of a power series development of a power series.
					template <class Series, class ArgsTuple>
					static size_t power_series_iterations(const Series &, const int &, const int &, const ArgsTuple &) {
						piranha_throw(value_error,"null truncator cannot provide number of iterations for power series");
					}
					bool is_effective() const {
						return false;
					}
			};
	};
}

#endif
