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

#ifndef PIRANHA_Q_POWER_CACHE_H
#define PIRANHA_Q_POWER_CACHE_H

#include <boost/unordered_map.hpp>
#include <vector>

#include "exceptions.h"
#include "mp.h"
#include "settings.h"

namespace piranha
{
	// Computed values are stored in a hash
	// map and re-used to calculate new values. To use it, simply construct and request the
	// value with operator[].
	template <class T, class ArithmeticFunctor>
	class q_power_cache
	{
		protected:
			typedef boost::unordered_map<mp_rational,T> container_type;
			typedef typename container_type::iterator iterator;
			typedef ArithmeticFunctor arith_functor_type;
			q_power_cache():m_container(),m_arith_functor() {}
		public:
			q_power_cache(const T &x): m_container(), m_arith_functor()
			{
				init(x);
			}
			const T &operator[](const mp_rational &q)
			{
				iterator it = m_container.find(q);
				if (it == m_container.end()) {
					if (q == 0) {
						m_container[mp_rational(0)] = m_arith_functor.pow(T(),0);
						return m_container[mp_rational(0)];
					} else {
						return insert_new(q);
					}
				} else {
					return it->second;
				}
			}
		private:
			void init(const T &x)
			{
				m_container[mp_rational(1)] = x;
			}
			T &insert_new(const mp_rational &q)
			{
				m_container[q] = m_arith_functor.pow(m_container[mp_rational(1)],q);
				return m_container[q];
			}
		protected:
			container_type			m_container;
			const arith_functor_type	m_arith_functor;
	};
}

#endif
