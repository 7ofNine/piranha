/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *   Copyright (C) 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.
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

// NOTE: this allocator is modified from the "malloc_allocator" memory allocator
// from GCC 4.2.0. Original copyright notice follows.

// ORIGINAL COPYRIGHT BEGINS --------------------------------------------------
// Copyright (C) 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
// ORIGINAL COPYRIGHT ENDS ----------------------------------------------------

#ifndef PIRANHA_MM_MALLOC_ALLOCATOR_H
#define PIRANHA_MM_MALLOC_ALLOCATOR_H

#include <cstdlib>
#include <new>
#include <bits/functexcept.h>
#include <mm_malloc.h>

namespace piranha
{

  using std::size_t;
  using std::ptrdiff_t;

/// An allocator that uses _mm_malloc.
/**
 * All allocation calls _mm_malloc, all deallocation calls _mm_free.
 */
  template<typename _Tp, size_t _Alignment>
    class mm_malloc_allocator
  {
    public:
      typedef size_t     size_type;
      typedef ptrdiff_t  difference_type;
      typedef _Tp*       pointer;
      typedef const _Tp* const_pointer;
      typedef _Tp&       reference;
      typedef const _Tp& const_reference;
      typedef _Tp        value_type;

      template<typename _Tp1>
        struct rebind
        { typedef mm_malloc_allocator<_Tp1,_Alignment> other; };

      mm_malloc_allocator() throw() { }

      mm_malloc_allocator(const mm_malloc_allocator&) throw() { }

      template<typename _Tp1>
        mm_malloc_allocator(const mm_malloc_allocator<_Tp1,_Alignment>&) throw() { }

      ~mm_malloc_allocator() throw() { }

      pointer
        address(reference __x) const { return &__x; }

      const_pointer
        address(const_reference __x) const { return &__x; }

// NB: __n is permitted to be 0.  The C++ standard says nothing
// about what the return value is when __n == 0.
      pointer
        allocate(size_type __n, const void* = 0)
      {
        std::cout << "Allocating\n";
        if (__builtin_expect(__n > this->max_size(), false))
          std::__throw_bad_alloc();

        pointer __ret = static_cast<_Tp*>(_mm_malloc(__n * sizeof(_Tp),_Alignment));
        if (!__ret)
          std::__throw_bad_alloc();
        return __ret;
      }

// __p is not permitted to be a null pointer.
      void
        deallocate(pointer __p, size_type)
        { _mm_free(static_cast<void*>(__p)); }

      size_type
        max_size() const throw()
        { return size_t(-1) / sizeof(_Tp); }

// _GLIBCXX_RESOLVE_LIB_DEFECTS
// 402. wrong new expression in [some_] allocator::construct
      void
        construct(pointer __p, const _Tp& __val)
        { ::new(__p) value_type(__val); }

      void
        destroy(pointer __p) { __p->~_Tp(); }
  };

  template<typename _Tp,size_t _Alignment>
    inline bool
    operator==(const mm_malloc_allocator<_Tp,_Alignment>&, const mm_malloc_allocator<_Tp,_Alignment>&)
    { return true; }

  template<typename _Tp,size_t _Alignment>
    inline bool
    operator!=(const mm_malloc_allocator<_Tp,_Alignment>&, const mm_malloc_allocator<_Tp,_Alignment>&)
    { return false; }
}
#endif
