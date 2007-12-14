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

#ifndef PIRANHA_CFS_H
#define PIRANHA_CFS_H

#include "../bits/poisson_series/coefficients/double_cf.h"
#include "../bits/poisson_series/norm_index.h"
#include "../bits/poisson_series/generic_fs.h"
#include "../bits/poisson_series/terms/simple_term.h"
#include "../bits/poisson_series/trigonometric_parts/trig_sparse_array.h"

namespace piranha
{
  typedef generic_fs<double_cf,trig_sparse_array,simple_term,norm_index> cfs;
}
#endif
