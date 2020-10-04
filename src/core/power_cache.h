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


#include <vector>
#include <unordered_map>

#include <boost/unordered_map.hpp>

#include "exceptions.h"
#include "settings.h"

namespace piranha
{
	// Simple power cache: will calculate power of objects and store them internally for reuse.
	template <class Argument, class T, class ArithmeticFunctor>
	class PowerCache
	{
		protected:

			typedef boost::unordered_map<T, Argument> ContainerType;
			typedef typename ContainerType::iterator  Iterator;

			PowerCache():container(), arithmeticFunctor() {}

		public:

			PowerCache(const Argument &x): container(), arithmeticFunctor()
			{
				container[T(1)] = x;
			}


			PowerCache(const Argument &x_1, const Argument &inv_x): container(), arithmeticFunctor()
			{
				container[T(1)]  = x_1;
				container[T(-1)] = inv_x;
			}


			const Argument &operator[](const T &x)
			{
				Iterator it = container.find(x);
				if (it == container.end()) 
				{
					if (x == T(0)) 
					{
						container[T(0)] = arithmeticFunctor.pow(Argument(), 0);
						return container[T(0)];

					} else 
					{
						return insertNew(x);
					}
				} else 
				{
					return it->second;
				}
			}


		private:

			Argument &insertNew(const T &x)
			{
				if (x < T(0)) 
				{
					Iterator it = container.find(T(-1));
					if (it != container.end()) 
					{
						container[x] = arithmeticFunctor.pow(container[T(-1)], -x);
						return container[x];
					}
				}

				container[x] = arithmeticFunctor.pow(container[T(1)], x);
				return container[x];
			}

		protected:
			ContainerType		    container;
			const ArithmeticFunctor	arithmeticFunctor;
	};
}

#endif
