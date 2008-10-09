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

#ifndef PIRANHA_PLAIN_SERIES_MULTIPLIER_H
#define PIRANHA_PLAIN_SERIES_MULTIPLIER_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/type_traits/is_same.hpp> // For key type detection.
#include <vector>

#include "../p_assert.h"
#include "../settings.h"
#include "base_series_multiplier.h"

namespace piranha
{
	/// Plain series multiplier.
	template <class Series1, class Series2, class ArgsTuple, template <class> class Truncator>
	class plain_series_multiplier: base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
				plain_series_multiplier<Series1, Series2, ArgsTuple, Truncator> >
	{
			typedef base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
			plain_series_multiplier<Series1, Series2, ArgsTuple, Truncator> > ancestor;
		public:
			// These typedefs are public because truncators may want to use them.
			/// Alias for term type of first input series and return value series.
			typedef typename Series1::term_type term_type1;
			/// Alias for term type of second input series.
			typedef typename Series2::term_type term_type2;
			/// Alias for the truncator type.
			typedef Truncator<plain_series_multiplier> truncator_type;
			plain_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
					ancestor::base_series_multiplier(s1, s1, retval, args_tuple) {}
			/// Perform multiplication.
			/**
			 * Method called by piranha::series_multiplication. Internally it calls perform_plain_multiplication.
			 */
			void perform_multiplication() {
				ancestor::perform_plain_multiplication();
			}
	};
}

#endif
