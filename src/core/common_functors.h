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

#include "exceptions.h"

// NOTE: maybe better move these and the caches into own subdirectory.

namespace piranha
{
	// TODO: use std::pow instead of pow method and inv().
	template <class T>
	struct named_series_arithmetics
	{
		T inv(const T &orig) const
		{
			return orig.pow(-1);
		}
		void multiply(T &orig, const T &other) const
		{
			orig *= other;
		}
		template <class U>
		T pow(const T &orig, const U &y) const
		{
			return orig.pow(y);
		}
	};

	template <class T, class ArgsTuple>
	struct base_series_arithmetics {
		base_series_arithmetics():m_argsTuple(0) {}
		T inv(const T &orig) const
		{
			piranha_assert(m_argsTuple);
			return orig.base_pow(-1,*m_argsTuple);
		}
		void multiply(T &orig, const T &other) const
		{
			piranha_assert(m_argsTuple);
			orig.base_mult_by(other,*m_argsTuple);
		}
		template <class U>
		T pow(const T &orig, const U &y) const
		{
			piranha_assert(m_argsTuple);
			return orig.base_pow(y,*m_argsTuple);
		}
		mutable ArgsTuple const *m_argsTuple;
	};

	struct ei_sub_functor {
		template <class RetSeries, class Element, class PosTuple, class SubCaches,
			class ArgsTuple>
		static RetSeries run(const Element &e, const PosTuple &pos_tuple,
			SubCaches &sub_caches, const ArgsTuple &argsTuple) {
			return e.template ei_sub<RetSeries>(pos_tuple, sub_caches, argsTuple);
		}
	};
}

#endif
