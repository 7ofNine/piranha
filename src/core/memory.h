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

#include <cstdlib> // For malloc.
#include <cstring> // For memcpy.
#include <exception> // For standard bad_alloc exception.

#include "config.h" // For unlikely(), aligned malloc, etc.
#include "math.h" // For lg to detect that memory alignment is a power of 2.

namespace piranha
{
	/// Low level memory allocation function.
	/**
	 * Thin wrapper around malloc(), will throw an instance of std::bad_alloc if allocation fails.
	 */
	inline void *piranha_malloc(const size_t &size)
	{
		void *retval = malloc(size);
		if (unlikely(retval == NULL)) {
			throw std::bad_alloc();
		}
		return retval;
	}

	template <int Alignment>
	inline void memory_alignment_checks() {
		p_static_check(Alignment > 0, "Memory alignment must be strictly positive.");
		// Test that Alignment is a multiple of the size of pointers.
		p_static_check(Alignment % sizeof(void *) == 0, "Memory alignment must be a multiple of sizeof(void *).");
		// Test that Alignment is at least as big as the size of pointers.
		p_static_check(Alignment >= sizeof(void *), "Memory alignment must be equal to or greater than sizeof(void *).");
		// Test that Alignment is a power of 2.
		p_static_check(lg<Alignment>::value > 0, "Memory alignment must be a multiple of 2.");
	}

	/// Low level memory allocation function supporting alignment specification.
	/**
	 * Thin wrapper around malloc(), will throw an instance of std::bad_alloc if allocation fails.
	 */
	template <int Alignment>
	inline void *piranha_malloc(const size_t &size)
	{
		memory_alignment_checks<Alignment>();
		void *ptr;
#ifdef _PIRANHA_WIN32
		ptr = _aligned_malloc(size,Alignment);
		if (unlikely(ptr == NULL)) {
#else
		if (unlikely(posix_memalign(&ptr,Alignment,size) != 0)) {
#endif
			throw std::bad_alloc();
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

	// To be used in conjunction with the aligning piranha_malloc.
	template <int Alignment>
	inline void piranha_free(void *ptr) {
		memory_alignment_checks<Alignment>();
#ifdef _PIRANHA_WIN32
		_aligned_free(ptr);
#else
		free(ptr);
#endif
	}
}

#endif
