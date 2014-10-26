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

#ifndef PIRANHA_POWER_CACHE_H
#define PIRANHA_POWER_CACHE_H

#include <boost/unordered_map.hpp>
#include <vector>

#include "exceptions.h"
#include "settings.h"

namespace piranha
{
	// Simple power cache: will calculate power of objects and store them internally for reuse.
	template <class Argument, class T, class ArithmeticFunctor>
	class PowerCache
	{
		protected:

			typedef boost::unordered_map<T, Argument> container_type;
			typedef typename container_type::iterator iterator;

			PowerCache():m_container(), arithmeticFunctor() {}

		public:

			PowerCache(const Argument &x): m_container(), arithmeticFunctor()
			{
				m_container[T(1)] = x;
			}


			PowerCache(const Argument &x_1, const Argument &inv_x): m_container(), arithmeticFunctor()
			{
				m_container[T(1)]  = x_1;
				m_container[T(-1)] = inv_x;
			}


			const Argument &operator[](const T &x)
			{
				iterator it = m_container.find(x);
				if (it == m_container.end()) 
				{
					if (x == 0) 
					{
						m_container[T(0)] = arithmeticFunctor.pow(Argument(),0);
						return m_container[T(0)];

					} else 
					{
						return insert_new(x);
					}
				} else 
				{
					return it->second;
				}
			}


		private:

			Argument &insert_new(const T &x)
			{
				if (x < 0) 
				{
					iterator it = m_container.find(T(-1));
					if (it != m_container.end()) 
					{
						m_container[x] = arithmeticFunctor.pow(m_container[T(-1)],-x);
						return m_container[x];
					}
				}

				m_container[x] = arithmeticFunctor.pow(m_container[T(1)],x);
				return m_container[x];
			}

		protected:
			container_type		    m_container;
			const ArithmeticFunctor	arithmeticFunctor;
	};
}

#endif
