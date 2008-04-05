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

#include <cstring>                                // For memcpy.
#include <ext/pool_allocator.h>
#include <gmp.h>
#include <gmpxx.h>

namespace piranha
{
  struct memory
  {
    static __gnu_cxx::__pool_alloc<char>  pool_allocator;
  };
}


inline void *mp_alloc(size_t size)
{
  return static_cast<void *>(piranha::memory::pool_allocator.allocate(size));
}


inline void mp_free(void *ptr, size_t size)
{
  piranha::memory::pool_allocator.deallocate((char *)ptr,size);
}


inline void *mp_realloc(void *ptr, size_t old_size, size_t new_size)
{
  void *retval=mp_alloc(new_size);
  memcpy(retval,ptr,old_size);
  mp_free(ptr,old_size);
  return retval;
}
#endif
