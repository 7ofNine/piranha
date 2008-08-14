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

#ifndef PIRANHA_ATOMIC_COUNTER_GENERIC_H
#define PIRANHA_ATOMIC_COUNTER_GENERIC_H

#ifdef _PIRANHA_MT
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#endif

#include "base_classes/base_atomic_counter.h"

namespace piranha
{
	template <class IntType>
	class atomic_counter_generic: public base_atomic_counter<IntType,atomic_counter_generic<IntType> >
	{
			typedef base_atomic_counter<IntType,atomic_counter_generic<IntType> > ancestor;
		public:
			atomic_counter_generic():ancestor::base_atomic_counter()
#ifdef _PIRANHA_MT
				,m_mutex()
#endif
			{}
			template <class IntType2>
			atomic_counter_generic &operator+=(const IntType2 &n) {
#ifdef _PIRANHA_MT
				boost::lock_guard<boost::mutex> lock(m_mutex);
#endif
				this->m_value += n;
				return *this;
			}
			template <class IntType2>
			atomic_counter_generic &operator-=(const IntType2 &n) {
#ifdef _PIRANHA_MT
				boost::lock_guard<boost::mutex> lock(m_mutex);
#endif
				this->m_value -= n;
				return *this;
			}
#ifdef _PIRANHA_MT
		private:
			boost::mutex m_mutex;
#endif
	};
}

#endif
