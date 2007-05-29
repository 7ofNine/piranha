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

#ifndef PIRANHA_GSP_H
#define PIRANHA_GSP_H

#include "double_cf.h"
#include "polynomial.h"
#include "ps.h"
#include "ps_term.h"
#include "trig_array.h"

namespace piranha
  {
  //typedef ps_term<polynomial<monomial_gmp_array<double_cf> >,trig_array> gsp_term;
  //typedef ps_term<complex_version<double_cf>,trig_array> npc_term;
  typedef ps<polynomial<double_cf>,trig_array,default_ps_index> gsp;
  //typedef ps<complex_version<double_cf>,trig_array> npc;
}

#endif
