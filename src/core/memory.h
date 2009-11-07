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
#include <cstddef>
#include <stdexcept>
#include <memory>

#include "atomic_counters/atomic_counters.h" // For counting allocator.
#include "base_classes/base_counting_allocator.h"
#include "config.h" // For unlikely().
#include "exceptions.h"
#include "settings.h"

namespace piranha
{
	/// STL-compatible allocator that decorates an existing allocator by adding a counting mechanism for the number of allocated bytes.
	/**
	 * The counting mechanism is thread-safe through the use of atomic operations.
	 */
	template <class T, class Allocator>
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
				const std::size_t add = n * sizeof(value_type), cur = m_counter.get_value(),
					l = settings::get_memory_limit();
				// Formulate in this way in order to avoid bogus values when doing l - add
				// (which is unsigned arithmetic).
				if (unlikely(add > l || cur > l - add)) {
					piranha_throw(memory_error,"memory limit reached");
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
				return boost::integer_traits<size_type>::const_max / sizeof(value_type);
			}
			void construct(pointer p, const value_type &val)
			{
				m_alloc.construct(p,val);
			}
			void destroy(pointer p)
			{
				m_alloc.destroy(p);
			}
		private:
			alloc m_alloc;
	};

	template<class T, class Allocator>
	inline bool operator==(const counting_allocator<T,Allocator> &, const counting_allocator<T,Allocator> &)
	{
		return true;
	}

	template<class T, class Allocator>
	inline bool operator!=(const counting_allocator<T,Allocator> &, const counting_allocator<T,Allocator> &)
	{
		return false;
	}

	template <class T>
	class std_counting_allocator: public counting_allocator<T,std::allocator<char> > {};
}

#endif
