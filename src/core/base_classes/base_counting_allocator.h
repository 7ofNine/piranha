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

#ifndef PIRANHA_BASE_COUNTING_ALLOCATOR_H
#define PIRANHA_BASE_COUNTING_ALLOCATOR_H

#include "../config.h"

#include <atomic>

namespace piranha
{
	class PIRANHA_VISIBLE BaseCountingAllocator   //this way all classes derived from this one share the counter
	{
	public:
		using CounterType = std::atomic_size_t;
	
	
		static CounterType::value_type count() noexcept
		{
			return counter;
		}
	protected:

			static inline CounterType counter = 0;
	};
}

#endif
