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

#ifndef PIRANHA_STATS_H
#define PIRANHA_STATS_H

#include "poisson_series/base_classes/base_pseries_def.h"

namespace piranha
{
/// Piranha-specific statistics class.
  class stats
  {
    template <__PIRANHA_BASE_PS_TP_DECL>
      friend class base_pseries;
    public:
      static double pack_ratio();
    private:
      static void insert()
      {
        ++total_insertions_;
      }
      static void pack()
      {
        insert();
        ++packed_insertions_;
      }
    private:
      static double  total_insertions_;
      static double  packed_insertions_;
  };
}
#endif                                            // PIRANHA_STATS_H
