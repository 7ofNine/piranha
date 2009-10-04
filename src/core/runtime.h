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

#ifndef PIRANHA_RUNTIME_H
#define PIRANHA_RUNTIME_H

#include <cstddef>

#include "atomic_counters/atomic_counters.h"

namespace piranha
{
	// Class that holds runtime info.
	class runtime {
		public:
			class register_threads {
				public:
					register_threads(const std::size_t &n):m_n(n)
					{
						m_n_current_threads += n;
					}
					~register_threads()
					{
						m_n_current_threads -= m_n;
					}
				private:
					const std::size_t m_n;
			};
			static std::size_t get_n_cur_threads()
			{
				return m_n_current_threads.get_value();
			}
		private:
			static atomic_counter_size_t m_n_current_threads;
	};
}

#endif
