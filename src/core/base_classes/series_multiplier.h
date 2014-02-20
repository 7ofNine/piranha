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

#ifndef PIRANHA_SERIES_MULTIPLIER_H
#define PIRANHA_SERIES_MULTIPLIER_H

#include "../base_classes/base_series_multiplier.h"

namespace piranha
{
	/// Generic series multiplier.
	class series_multiplier
	{
		public:
			template <class Series1, class Series2, class ArgsTuple, class Truncator>
			class get_type:
				public base_series_multiplier<Series1, Series2, ArgsTuple, Truncator, get_type<Series1, Series2, ArgsTuple, Truncator> >
			{
					typedef base_series_multiplier<Series1, Series2, ArgsTuple, Truncator,
						get_type<Series1, Series2, ArgsTuple, Truncator> > ancestor;
				public:
					// TODO: which of these are needed? Cleanup also elsewhere when parallel is implemented.
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple ArgsTupleType;
					typedef typename Truncator::template get_type<Series1, Series2, ArgsTuple> truncator_type;

					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &argsTuple):ancestor(s1, s2, retval, argsTuple) {}

					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication()
					{
						// Cache term pointers.
						this->cache_terms_pointers(this->m_s1, this->m_s2);
						this->perform_plain_multiplication();
					}
			};
	};
}

#endif
