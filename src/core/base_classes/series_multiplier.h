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
#include "../base_classes/null_truncator.h"

namespace piranha
{
	/// Generic series multiplier.
	// NOTE: isn't truncator's "accept" used here?
	class series_multiplier
	{
		public:
			template <class Series1, class Series2, class ArgsTuple, class Truncator>
			class get_type:
				public base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
				get_type<Series1, Series2, ArgsTuple, Truncator> >
			{
					typedef base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
						get_type<Series1, Series2, ArgsTuple, Truncator> > ancestor;
// 					typedef typename Series1::const_iterator const_iterator1;
// 					typedef typename Series2::const_iterator const_iterator2;
// 					typedef typename ancestor::term_type1 term_type1;
// 					typedef typename ancestor::term_type2 term_type2;
// 					typedef typename term_type1::cf_type cf_type1;
// 					typedef typename term_type2::cf_type cf_type2;
// 					typedef typename term_type1::key_type key_type;
				public:
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<get_type> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):ancestor(s1, s2, retval, args_tuple)
					{}
					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication()
					{
						// Build the truncator.
						const truncator_type trunc(*this);
						// Use the selected truncator only if it really truncates, otherwise use the
						// null truncator.
						if (trunc.is_effective()) {
							ll_perform_multiplication(trunc);
						} else {
							ll_perform_multiplication(null_truncator::template get_type<get_type>(*this));
						}
					}
				private:
					template <class GenericTruncator>
					void ll_perform_multiplication(const GenericTruncator &trunc)
					{
						this->perform_plain_multiplication(trunc);
					}
			};
	};
}

#endif
