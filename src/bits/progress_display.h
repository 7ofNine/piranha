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

#ifndef PIRANHA_PROGRESS_DISPLAY_H
#define PIRANHA_PROGRESS_DISPLAY_H

#include <boost/progress.hpp>
#include <boost/scoped_ptr.hpp>

#include "config.h" // For _PIRANHA_DISPLAY_PROGRESS_MAX_N.

namespace piranha
{
  template <bool Active>
    class progress_display
  {
    public:
      progress_display(const size_t &n):active(n>_PIRANHA_DISPLAY_PROGRESS_MAX_N),pd(0)
      {
        switch (active)
        {
          case true:
            pd.reset(new boost::progress_display(n));
          case false:
            ;
        }
      }
      void operator++()
      {
        switch (active)
        {
          case true:
            ++(*pd);
          case false:
            ;
        }
      }
      unsigned long int count() const
      {
        switch (active)
        {
          case true:
            return pd->count();
          default:
            return 0;
        }
      }
    private:
      const bool                                  active;
      boost::scoped_ptr<boost::progress_display>  pd;
  };

  template <>
    class progress_display<false>
  {
    public:
      progress_display(const size_t &) {}
      void operator++() {}
      unsigned long int count() const
      {
        return 0;
      }
  };
}

#endif
