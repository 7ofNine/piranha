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

#ifndef PIRANHA_SNP_H
#define PIRANHA_SNP_H

#ifdef _PIRANHA_SSE2

#include "double_cf.h"
#include "pool_allocator.h"
#include "norm_index.h"
#include "ps.h"
#include "simple_term.h"
#include "trig_simd_array.h"

namespace piranha
{
  template <int N>
    struct snp
  {
    typedef ps<double_cf,trig_simd_array<N>,simple_term,norm_based_index,pool_allocator<char,16> > type;
  };
}

#endif

#endif