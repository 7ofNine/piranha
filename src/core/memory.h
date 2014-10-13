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

#include <boost/cstdint.hpp> // For uint8_t.
#include <boost/integer_traits.hpp> // For max allocatable number of objects.
#include <boost/type_traits/is_same.hpp> // For type mismatch identification in the counting allocator.
#include <stdint.h>
#include <cstdlib>
#include <cstddef>
#include <memory>
#include <new> // For std::bad_alloc.

#include "atomic_counters/atomic_counters.h" // For counting allocator.
#include "base_classes/base_counting_allocator.h"
#include "config.h" // For unlikely().
#include "exceptions.h"
#include "integer_typedefs.h"
#include "settings.h"

namespace piranha
{
	// TODO: get rid of the T parameter, std_counting_allocator and friends.
	/// STL-compatible allocator that decorates an existing allocator by adding a counting mechanism for the number of allocated bytes.
	/**
	 * The counting mechanism is thread-safe through the use of atomic operations. The counting is approximate, in the sense that while concurrent
	 * allocations are taking place it is not guaranteed that in each moment the count reflects exactly the number of allocated bytes.
	 */
	template <class T, class Allocator>
	class CountingAllocator: public base_counting_allocator
	{
			template <class U, class Allocator2>
			friend class CountingAllocator;
			typedef typename Allocator::template rebind<T>::other alloc;
			PIRANHA_STATIC_CHECK((boost::is_same<T,typename alloc::value_type>::value), "Type mismatch in counting allocator.");
		public:
			typedef typename alloc::size_type size_type;
			typedef typename alloc::difference_type difference_type;
			typedef typename alloc::pointer pointer;
			typedef typename alloc::const_pointer const_pointer;
			typedef typename alloc::reference reference;
			typedef typename alloc::const_reference const_reference;
			typedef typename alloc::value_type value_type;
			template <class U>
			struct rebind {
				typedef CountingAllocator<U,Allocator> other;
			};
			CountingAllocator():m_alloc() {}
			CountingAllocator(const CountingAllocator &a):m_alloc(a.m_alloc) {}
			template <class U>
			CountingAllocator(const CountingAllocator<U,Allocator> &a):m_alloc(a.m_alloc) {}
			pointer address(reference x)
			{
				return m_alloc.address(x);
			}
			const_pointer address(const_reference x) const
			{
				return m_alloc.address(x);
			}
			pointer allocate(const size_type &n, const void *hint = 0)
			{
				// TODO: guard overflow here?
				const std::size_t add = n * sizeof(value_type), cur = m_counter.get_value(),
					l = settings::get_memory_limit();
				// Formulate in this way in order to avoid bogus values when doing l - add
				// (which is unsigned arithmetic).
				if (unlikely(add > l || cur > l - add)) {
					PIRANHA_THROW(memory_error,"memory limit reached");
				}
				pointer retval = m_alloc.allocate(n,hint);
				m_counter += add;
				return retval;
			}
			void deallocate(pointer p, const size_type &n)
			{
				m_alloc.deallocate(p,n);
				m_counter -= n * sizeof(value_type);
			}
			size_type max_size() const
			{
				return m_alloc.max_size();
			}
			void construct(pointer p, const value_type &val)
			{
				m_alloc.construct(p,val);
			}
			void destroy(pointer p)
			{
				m_alloc.destroy(p);
			}
			bool operator==(const CountingAllocator &c) const
			{
				return (m_alloc == c.m_alloc);
			}
			bool operator!=(const CountingAllocator &c) const
			{
				return (m_alloc != c.m_alloc);
			}
		private:
			alloc m_alloc;
	};

	template <class T>
	class std_counting_allocator: public CountingAllocator<T,std::allocator<char> > {};

	// The following allocator was slightly adapted from the XVID project. Original copyright notice follows.
	/*****************************************************************************
	*
	*  XVID MPEG-4 VIDEO CODEC
	*  - Aligned Memory Allocator -
	*
	*  Copyright(C) 2002-2003 Edouard Gomez <ed.gomez@free.fr>
	*
	*  This program is free software ; you can redistribute it and/or modify
	*  it under the terms of the GNU General Public License as published by
	*  the Free Software Foundation ; either version 2 of the License, or
	*  (at your option) any later version.
	*
	*  This program is distributed in the hope that it will be useful,
	*  but WITHOUT ANY WARRANTY ; without even the implied warranty of
	*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	*  GNU General Public License for more details.
	*
	*  You should have received a copy of the GNU General Public License
	*  along with this program ; if not, write to the Free Software
	*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
	*
	* $Id: mem_align.c,v 1.16 2004/03/22 22:36:24 edgomez Exp $
	*
	****************************************************************************/
	/// STL-compatible allocator that aligns memory to a specific boundary N.
	/**
	 * N must be in the [0,256[ range. If N == 0, the alignment will be undefined.
	 * Implementation uses std::malloc(). Slightly adapted from http://cvs.xvid.org/cvs/viewvc.cgi/xvidcore/src/utils/mem_align.c?view=log.
	 */
	template <class T, int Alignment>
	class align_mallocator {

			PIRANHA_STATIC_CHECK(Alignment >= 0 && Alignment < 256,"Invalid alignment value requested.");

		public:

			typedef std::size_t	size_type;
			typedef std::ptrdiff_t	difference_type;
			typedef T *		pointer;
			typedef const T*	const_pointer;
			typedef T&		reference;
			typedef const T&	const_reference;
			typedef T		value_type;
			template <class U>
			struct rebind {
				typedef align_mallocator<U,Alignment> other;
			};

			align_mallocator() {}

			align_mallocator(const align_mallocator&) {}

			template <class U>
			align_mallocator(const align_mallocator<U,Alignment> &) {}

			~align_mallocator() {}

			pointer address(reference x) const
			{
				return &x;
			}

			const_pointer address(const_reference x) const
			{
				return &x;
			}

			pointer allocate(const size_type &n, const void * = 0)
			{
				// Just return 0 if no space is required.
				if (unlikely(!n)) 
                {
					return pointer(0);
				}
				if (unlikely(n > max_size())) 
                {
					throw std::bad_alloc();
				}
				boost::uint8_t *mem_ptr;
				if (!Alignment) 
                {
					// We have not to satisfy any alignment.
					if ((mem_ptr = (boost::uint8_t *)std::malloc(n * sizeof(value_type) + 1)) != 0) 
                    {
						// Store (mem_ptr - "real allocated memory") in *(mem_ptr-1).
						*mem_ptr = (boost::uint8_t)1;
						// Return the mem_ptr pointer.
						return ((pointer)(mem_ptr + 1));
					}
				} else 
                {
					boost::uint8_t *tmp;
					// Allocate the required size memory + alignment so we
					// can realign the data if necessary.
					if ((tmp = (boost::uint8_t *)std::malloc(n * sizeof(value_type) + Alignment)) != 0) 
                    {
						// Align the tmp pointer.
						mem_ptr = (boost::uint8_t *)((ptr_uint_t)(tmp + Alignment - 1) & (~(ptr_uint_t)(Alignment - 1)));
						// Special case where malloc has already satisfied the alignment
						// We must add alignment to mem_ptr because we must store
						// (mem_ptr - tmp) in *(mem_ptr-1)
						// If we do not add alignment to mem_ptr then *(mem_ptr-1) points
						// to a forbidden memory space.
						if (mem_ptr == tmp) 
                        {
							mem_ptr += Alignment;
						}
						// (mem_ptr - tmp) is stored in *(mem_ptr-1) so we are able to retrieve
						// the real malloc block allocated and free it in deallocation.
						*(mem_ptr - 1) = (boost::uint8_t)(mem_ptr - tmp);
						// Return the aligned pointer
						return ((pointer)mem_ptr);
					}
				}

				throw std::bad_alloc();
			}

			void deallocate(pointer p, const size_type &)
			{
				// Take care of zero allocation.
				if (unlikely(!p)) 
                {
					return;
				}
				// Aligned pointer.
				boost::uint8_t *ptr = (boost::uint8_t *)p;
				// *(ptr - 1) holds the offset to the real allocated block
				// we sub that offset os we free the real pointer.
				ptr -= *(ptr - 1);
				// Free the memory.
				std::free(ptr);
			}

			size_type max_size() const
			{
				// TODO: replace with boost numerical limits for integers.
				return std::size_t(-1) / sizeof(value_type);
			}

			void construct(pointer p, const T &value)
			{
				::new((void *)p) value_type(value);
			}

			void destroy(pointer p)
			{
				p->~T();
			}
	};

	template <class T, int Alignment>
	inline bool operator==(const align_mallocator<T, Alignment> &, const align_mallocator<T, Alignment> &)
	{
		return true;
	}

	template <class T, int Alignment>
	inline bool operator!=(const align_mallocator<T, Alignment> &, const align_mallocator<T, Alignment> &)
	{
		return false;
	}
}

#endif
