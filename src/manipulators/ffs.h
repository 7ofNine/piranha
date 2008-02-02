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

#ifndef PIRANHA_FFS_H
#define PIRANHA_FFS_H

#include "../bits/poisson_series/coefficients/double_cf.h"
#include "../bits/poisson_series/norm_index.h"
#ifdef _PIRANHA_SSE2
#include "../bits/pool_allocator.h"
#endif
#include "../bits/poisson_series/generic_fs.h"
#include "../bits/poisson_series/terms/simple_term.h"
#include "../bits/poisson_series/trigonometric_parts/trig_fixed_array.h"

namespace piranha
{
  template <int N>
    struct ffs
  {
    typedef generic_fs<double_cf,
      trig_fixed_array<N,16,1>,
      simple_term,
      norm_index
#ifdef _PIRANHA_SSE2
      , pool_allocator<char,16>
#endif
    > type;
  };
}
#endif
