/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_BUFFER_H
#define PIRANHA_BUFFER_H

#include "piranha_malloc.h"

namespace piranha
{
  class buffer_init
  {
    public:
// Allocate 4MB for start.
// TODO: maybe specify this value somewhere?
#define _INIT_BUFFER_SIZE_MB 200
      buffer_init():size(_INIT_BUFFER_SIZE_MB*1024*1024),ptr(piranha_malloc(size))
      {
        std::cout << "Buffer set up, around " << size/(1024*1024) << " MBytes available." << std::endl;
      }
#undef _INIT_BUFFER_SIZE_MB
      ~buffer_init()
      {
        piranha_free(ptr);
      }
      void *pointer()
      {
        return ptr;
      }
/// Returns size in bytes.
      const size_t &g_size() const
      {
        return size;
      }
    private:
      size_t  size;
      void    *ptr;
  };

  class buffer
  {
    public:
      template <class T>
        static T *head()
      {
        return (T *)bi.pointer();
      }
      static const size_t &g_size()
      {
        return bi.g_size();
      }
      template <class T>
        static size_t n_elements()
      {
        return g_size()/sizeof(T);
      }
    private:
      static buffer_init  bi;
  };

  buffer_init buffer::bi;
}

#endif
