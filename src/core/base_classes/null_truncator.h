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
	class NullTruncator
	{
		public:
			template <class Series1, class Series2, class ArgsTuple>
			class GetType
			{
				public:
					typedef typename Series1::TermType TermType1;
					typedef typename Series2::TermType TermType2;
					typedef GetType Type;

					GetType(std::vector<TermType1 const *> &, std::vector<TermType2 const *> &, ArgsTuple const &) {}


					template <class Term1, class Term2>
					bool skip(Term1 const &, Term2 const &) const
					{
						return false;
					}


					// Limit of a power series development of a power series.
					template <class Series, class ArgsTuple2>
					static  int powerSeriesIterations(Series const &, int const, int const , ArgsTuple2 const &)
					{
						PIRANHA_THROW(value_error, "null truncator cannot provide number of iterations for power series.");
					}


					bool isEffective() const
					{
						return false;
					}
			};

			template <class Series, class ArgsTuple2>
			static std::vector<typename Series::TermType const *> getSortedPointerVector(Series const &, ArgsTuple2 const &)
			{
				PIRANHA_THROW(value_error,"null truncator cannot order series");
			}
	};
}

#endif
