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

#include "tass17.h"

namespace piranha
{
// Initialization of tass17's static variables
  bool tass17::loaded=false;
  bool tass17::has_deltas=false;

  const double tass17::m0_=3498.790;
  const double tass17::m6_=0.4225863977890E+04;

  lnp tass17::lambda4_;
  lnp tass17::lambda6_;

  lnp tass17::p6_;

  lnpc tass17::z6_;

  lnpc tass17::zeta6_;

  lnp tass17::dlambda1_;
  lnp tass17::dlambda2_;
  lnp tass17::dlambda3_;
  lnp tass17::dlambda4_;
  lnp tass17::dlambda5_;
  lnp tass17::dlambda6_;
  lnp tass17::dlambda8_;
}
