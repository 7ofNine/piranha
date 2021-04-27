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

#include "base_classes/base_counting_allocator.h"
#include "exceptions.h"
#include "integer_typedefs.h"
#include "settings.h"
#include "config.h"

#include <cstdlib>
#include <cstddef>
#include <memory>
#include <new> // For std::bad_alloc.
#include <type_traits>
#include <atomic>

namespace piranha
{
	// TODO: get rid of the T parameter, std_counting_allocator and friends.
	/// STL-compatible allocator that decorates an existing allocator by adding a counting mechanism for the number of allocated bytes.
	//
	// The counting mechanism is thread-safe through the use of atomic operations. The counting is approximate though, in the sense that while concurrent
	// allocations are taking place it is not guaranteed that at each moment the count read reflects exactly the number of allocated bytes.
	// It may change between reaeding count and actual allocation, i.e. the underlying Allocator could still throw an exception.
	// We could make it completely safe but for a performance penalty vs. a small risk of triggering the exception which could happen
	// anyways as long as we using an non pre-allocated heap or memory pool.
	//

	template <typename T, typename Allocator = std::allocator<T> >
	class CountingAllocator: public BaseCountingAllocator
	{
			template <class U, class Allocator2>
			friend class CountingAllocator;

			using AllocatorType = typename std::allocator_traits<Allocator>:: template rebind_alloc<T>; //rebind to template type
			using AllocInterface = std::allocator_traits<AllocatorType>;

			static_assert((std::is_same_v<T, typename AllocatorType::value_type>), "Type mismatch in counting allocator.");

		public:
			typedef typename AllocatorType::size_type size_type;
			typedef typename AllocatorType::difference_type difference_type;
			typedef typename AllocatorType::value_type value_type;
			
			CountingAllocator() = default;
			
			CountingAllocator(const CountingAllocator& a) = default;

			//CountingAllocator(CountingAllocator &&) = delete;

			// CountingAllocator& operator=(CountingAllocator &&) = delete;

			~CountingAllocator() = default;
			
			template <class U>
			CountingAllocator(const CountingAllocator<U, Allocator> &a):m_alloc(a.m_alloc) {}    // how does this work for different sized object being allocated??

			[[nodiscard]] constexpr T* allocate(const size_type n)
			{
				// TODO: guard overflow here?
				const std::size_t add      = n * sizeof(value_type);
                std::size_t       current = counter;                  // current allocated memory size
				std::size_t       memLimit = settings::get_memory_limit();
				// Formulate in this way in order to avoid bogus values when doing memLimit - add
				// (which is unsigned arithmetic).
				if (add > memLimit || current > memLimit - add) // reqquest bigger than memory limit or would  exceeded it
                {
					PIRANHA_THROW(memory_error, "memory limit reached");
				}

				T* retval = AllocInterface::allocate(m_alloc, n);
				counter += add;   // this is atomic
				return retval;
			}

			constexpr void deallocate(T* p, const size_type n)
			{
				AllocInterface::deallocate(m_alloc, p, n);
				counter -= n * sizeof(value_type);  //this is atomic
			}

			bool operator==(const CountingAllocator&) const = default;

	private:
			AllocatorType m_alloc;
	};

}

#endif
