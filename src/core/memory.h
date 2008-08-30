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

#include <boost/integer_traits.hpp> // For max allocatable number of objects.
#include <boost/type_traits/is_same.hpp> // For type mismatch identification in the counting allocator.
#include <cstdlib> // For malloc.
#include <cstring> // For memcpy.
#include <exception> // For standard bad_alloc exception.
#include <memory>

#include "atomic_counter.h" // For counting allocator.
#include "base_classes/base_counting_allocator.h"
#include "config.h" // For unlikely(), aligned malloc, visibility, etc.
#include "exceptions.h"
#include "math.h" // For lg to detect that memory alignment is a power of 2.
#include "portable_memalign.h"
#include "settings.h"

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
		// Test that Alignment is no greater than 128 (2**7). This is done beacuse portable memalign
		// uses a uint8 to store the alignment info.
		p_static_check(Alignment <= 128, "Memory alignment value is too big (>128)");
	}

	/// Low level memory allocation function supporting alignment specification.
	/**
	 * Thin wrapper around a platform-specific memory aligning function, will throw an
	 * instance of std::bad_alloc if allocation fails.
	 */
	template <int Alignment>
	inline void *piranha_malloc(const size_t &size)
	{
		memory_alignment_checks<Alignment>();
		void *ptr;
#if defined(_MSC_VER)
		ptr = _aligned_malloc(size,Alignment);
		if (unlikely(ptr == NULL)) {
#elif defined (_GNUC_) && ( !defined (_PIRANHA_MINGW) )
		if (unlikely(posix_memalign(&ptr,Alignment,size) != 0)) {
#else
		ptr = portable_memalign<Alignment>(size);
		if (unlikely(ptr == NULL)) {
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
#if defined(_MSC_VER)
		_aligned_free(ptr);
#elif defined (_GNUC_) && ( !defined (_PIRANHA_MINGW) )
		free(ptr);
#else
		portable_aligned_free(ptr);
#endif
	}

	/// An allocator that decorates another allocator by adding a counting mechanism for the number of allocated bytes.
	template<class T, class Allocator>
	class counting_allocator: public base_counting_allocator
	{
			typedef typename Allocator::template rebind<T>::other alloc;
			p_static_check((boost::is_same<T,typename alloc::value_type>::value), "Type mismatch in counting allocator.");
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
				typedef counting_allocator<U,Allocator> other;
			};
			counting_allocator():m_alloc() {}
			counting_allocator(const counting_allocator &):m_alloc() {}
			template <class U>
			counting_allocator(const counting_allocator<U,Allocator> &):m_alloc() {}
			~counting_allocator() {}
			pointer address(reference x) const {
				return m_alloc.address(x);
			}
			const_pointer address(const_reference x) const {
				return m_alloc.address(x);
			}
			pointer allocate(const size_type &n, const void *hint = 0) {
				const size_t add = n * sizeof(value_type), cur = static_cast<size_t>(m_counter), l = settings::memory_limit();
				// Cast to double so that we resolve the case in which cur+add overflows size_t.
				if (static_cast<double>(cur) + add > l) {
					throw out_of_memory();
				}
				pointer retval = m_alloc.allocate(n,hint);
				m_counter += add;
				return retval;
			}
			void deallocate(pointer p, const size_type &n) {
				m_alloc.deallocate(p,n);
				m_counter -= n * sizeof(value_type);
			}
			size_type max_size() const {
				return boost::integer_traits<size_type>::const_max/sizeof(value_type);
			}
			void construct(pointer p, const value_type &val) {
				m_alloc.construct(p,val);
			}
			void destroy(pointer p) {
				m_alloc.destroy(p);
			}
		private:
			alloc m_alloc;
	};

	template<class T, class Allocator>
	inline bool operator==(const counting_allocator<T,Allocator> &, const counting_allocator<T,Allocator> &) {
		return true;
	}

	template<class T, class Allocator>
	inline bool operator!=(const counting_allocator<T,Allocator> &, const counting_allocator<T,Allocator> &) {
		return false;
	}

	template <class T>
	class std_counting_allocator: public counting_allocator<T,std::allocator<char> > {};
}

#endif
