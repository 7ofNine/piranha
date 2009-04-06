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

#ifndef PIRANHA_COMMON_FUNCTORS_H
#define PIRANHA_COMMON_FUNCTORS_H

#include <algorithm>

#include "p_assert.h"

namespace piranha
{
	template <class T>
	struct member_hash_value
	{
		size_t operator()(const T &x) const {
			return x.hash_value();
		}
	};

	template <class T>
	struct member_swap
	{
		void operator()(T &a, T &b) const {
			a.swap(b);
		}
	};

	template <class T>
	struct std_swap
	{
		void operator()(T &a, T &b) const {
			std::swap(a,b);
		}
	};

	template <class T>
	struct named_series_arithmetics
	{
		T inv(const T &orig) const {
			return orig.inv();
		}
		void multiply(T &orig, const T &other) const {
			orig *= other;
		}
	};

	template <class T, class ArgsTuple>
	struct base_series_arithmetics
	{
		base_series_arithmetics():m_args_tuple(0)
		{}
		T inv(const T &orig) const {
			p_assert(m_args_tuple);
			return orig.base_inv(*m_args_tuple);
		}
		void multiply(T &orig, const T &other) const {
			p_assert(m_args_tuple);
			orig.base_mult_by(other,*m_args_tuple);
		}
		mutable ArgsTuple const *m_args_tuple;
	};

	struct ei_sub_functor {
		template <class RetSeries, class Element, class PosTuple, class SubCaches,
			class ArgsTuple>
		static RetSeries run(const Element &e, const PosTuple &pos_tuple,
			SubCaches &sub_caches, const ArgsTuple &args_tuple) {
			return e.template ei_sub<RetSeries>(pos_tuple, sub_caches, args_tuple);
		}
	};
}

#endif
