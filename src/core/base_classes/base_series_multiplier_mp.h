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

#ifndef PIRANHA_BASE_SERIES_MULTIPLIER_MP_H
#define PIRANHA_BASE_SERIES_MULTIPLIER_MP_H

#include <boost/tuple/tuple.hpp>

namespace piranha
{
	class base_insert_multiplication_result
	{
		protected:
			template <class SingleRes, class BaseMultiplier, class Truncator>
			static void insert_single_res(const SingleRes &res, BaseMultiplier &m, const Truncator &trunc) {
				const typename BaseMultiplier::mult_set::const_iterator it = m.m_set.find(res);
				switch (trunc.accept(res)) {
				case true:
					switch (it == m.m_set.end()) {
						// Not duplicate, insert it.
					case true:
						m.m_set.insert(res);
						break;
						// Duplicate, modify existing.
					case false:
						it->m_cf.add(res.m_cf, m.m_args_tuple);
					}
					break;
				case false:
					;
				}
			}
	};

	// Traverse the tuple of results of multiplication and insert each result into the multiplication result set.
	template <class ResultTuple>
	class insert_multiplication_result: public base_insert_multiplication_result
	{
		public:
			template <class BaseMultiplier, class Truncator>
			static void run(const ResultTuple &mult_res, BaseMultiplier &m, const Truncator &trunc) {
				insert_single_res(mult_res.get_head(), m, trunc);
				insert_multiplication_result<typename ResultTuple::tail_type>::run(mult_res.get_tail(), m, trunc);
			}
	};

	template <>
	class insert_multiplication_result<boost::tuples::null_type>
	{
		public:
			template <class BaseMultiplier, class Truncator>
			static void run(const boost::tuples::null_type &, BaseMultiplier &, const Truncator &) {}
	};
}

#endif
