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

#include "core/astro.h"

namespace piranha
{
  // Initialize values for constants.
  const double astro::J2000dot0_ = 2451545.0;
  const double astro::J1980dot0_ = 2444240.0;
  const double astro::JD_per_JY_ = 365.25;
  const double astro::seconds_per_JY_ = astro::JD_per_JY_*86400.0;
  const double astro::k_ = 0.01720209895;
  // From http://ssd.jpl.nasa.gov/?constants.
  const double astro::AU_ = 1.49597870691e11;
  // Source =
  // E. Myles Standish. "Report of the IAU WGAS Sub-group on Numerical Standards".
  // In Highlights of Astronomy, I.
  // Appenzeller, ed. Dordrecht: Kluwer Academic Publishers, 1995. (Complete report available
  // online: PostScript. Tables
  // from the report also available: Astrodynamic Constants and Parameters).
  const double astro::G_ = 6.67259e-11;
  const double astro::eps_0_ = 4.0909262968940374e-01;
}
