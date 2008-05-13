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

#ifndef PIRANHA_MEMORY_H
#define PIRANHA_MEMORY_H

#include <boost/static_assert.hpp>
#include <cstdlib> // For malloc.
#include <cstring> // For memcpy.
#include <exception> // For standard bad_alloc exception.
#include <gmp.h>
#include <gmpxx.h>

#include "config.h" // For unlikely(), aligned malloc, etc.
#include "math.h" // For lg to detect that memory alignment is a power of 2.

namespace piranha
{
	/// Low level memory allocation function.
	/**
	 * Thin wrapper around malloc(), will throw an instance of std::bad_alloc if allocation fails.
	 */
	inline void *piranha_malloc(const size_t &size) throw(std::bad_alloc)
	{
		void *retval = malloc(size);
		switch (unlikely(retval == NULL)) {
		case true:
			throw std::bad_alloc();
			break;
		case false:
			;
		}
		return retval;
	}

	/// Low level memory allocation function supporting alignment specification.
	/**
	 * Thin wrapper around malloc(), will throw an instance of std::bad_alloc if allocation fails.
	 */
	template <int Alignment>
	inline void *piranha_malloc(const size_t &size) throw(std::bad_alloc)
	{
		BOOST_STATIC_ASSERT(Alignment > 0);
		// Test that Alignment is a multiple of the size of pointers.
		BOOST_STATIC_ASSERT(Alignment % sizeof(void *) == 0);
		// Test that Alignment is at least as big as the size of pointers.
		BOOST_STATIC_ASSERT(Alignment >= sizeof(void *));
		// Test that Alignment is a power of 2.
		BOOST_STATIC_ASSERT(lg<Alignment>::value > 0);
		void *ptr;
		switch (unlikely(__ALIGNED_MALLOC(&ptr, Alignment, size) == 0)) {
		case true:
			throw std::bad_alloc();
			break;
		case false:
			;
		}
		return ptr;
	}

	/// Low level memory deallocation function.
	/**
	 * Thin wrapper around free().
	 */
	inline void piranha_free(void *ptr)
	{
		free(ptr);
	}
}

#endif
