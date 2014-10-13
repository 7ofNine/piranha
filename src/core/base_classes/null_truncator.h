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

#include <vector>

namespace piranha
{
	/// Truncator which does not truncate.
	class null_truncator
	{
		public:
			template <class Series1, class Series2, class ArgsTuple>
			class get_type
			{
				public:
					typedef typename Series1::term_type term_type1;
					typedef typename Series2::term_type term_type2;
					typedef get_type type;

					get_type(std::vector<term_type1 const *> &, std::vector<term_type2 const *> &, const ArgsTuple &) {}


					template <class Term1, class Term2>
					bool skip(const Term1 &, const Term2 &) const
					{
						return false;
					}


					// Limit of a power series development of a power series.
					template <class Series, class ArgsTuple2>
					static size_t powerSeriesIterations(const Series &, const int &, const int &, const ArgsTuple2 &)
					{
						PIRANHA_THROW(value_error,"null truncator cannot provide number of iterations for power series");
					}


					bool isEffective() const
					{
						return false;
					}
			};

			template <class Series, class ArgsTuple2>
			static std::vector<typename Series::term_type const *> getSortedPointerVector(const Series &, const ArgsTuple2 &)
			{
				PIRANHA_THROW(value_error,"null truncator cannot order series");
			}
	};
}

#endif
