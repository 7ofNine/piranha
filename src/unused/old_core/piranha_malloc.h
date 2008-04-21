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

#ifndef PIRANHA_PIRANHA_MALLOC_H
#define PIRANHA_PIRANHA_MALLOC_H

#include <cstdlib>                                // For malloc.

#ifdef _PIRANHA_SSE2
#ifdef _PIRANHA_PRIVATE_MM_MALLOC_H
#include "mm_malloc.h"
#else
#include <mm_malloc.h>
#endif
#endif

#include "config.h"                               // For unlikely().

namespace piranha
{
  /// Memory allocation function.
  /**
   * Wrapper around malloc() (or _mm_malloc() if SSE2 are used).
   */
  static __inline void *piranha_malloc(const size_t &size)
  {
#ifdef _PIRANHA_SSE2
    void *retval = _mm_malloc(size,16);
#else
    void *retval = malloc(size);
#endif
    if (unlikely(retval == NULL))
    {
      std::cout << "Error allocating memory, bailing out." << std::endl;
      std::abort();
    }
    return retval;
  }

  /// Memory deallocation function.
  /**
   * Wrapper around free() (or _mm_free() if SSE2 are used).
   */
  static __inline void piranha_free(void *ptr)
  {
#ifdef _PIRANHA_SSE2
    _mm_free(ptr);
#else
    free(ptr);
#endif
  }
}
#endif                                            // PIRANHA_PIRANHA_MALLOC_H
