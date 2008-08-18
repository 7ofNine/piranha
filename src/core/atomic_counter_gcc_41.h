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

#ifndef PIRANHA_ATOMIC_COUNTER_GCC_41_H
#define PIRANHA_ATOMIC_COUNTER_GCC_41_H

#include "base_classes/base_atomic_counter.h"

namespace piranha
{
	template <class IntType>
	class atomic_counter_gcc_41: public base_atomic_counter<IntType, atomic_counter_gcc_41<IntType> >
	{
			typedef base_atomic_counter<IntType, atomic_counter_gcc_41<IntType> > ancestor;
		public:
			atomic_counter_gcc_41():ancestor::base_atomic_counter() {}
			template <class IntType2>
			atomic_counter_gcc_41 &operator+=(const IntType2 &n) {
				__sync_add_and_fetch(&(this->m_value),static_cast<IntType>(n));
				return *this;
			}
			template <class IntType2>
			atomic_counter_gcc_41 &operator-=(const IntType2 &n) {
				__sync_sub_and_fetch(&(this->m_value),static_cast<IntType>(n));
				return *this;
			}
	};
}

#endif
