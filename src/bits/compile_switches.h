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

#ifndef PIRANHA_COMPILE_SWITCHES_H
#define PIRANHA_COMPILE_SWITCHES_H

namespace piranha
{
  /// Compile-time switches.
  /**
   * This struct reports some useful compile-time switches in form of boolean flags.
   */
  struct compile_switches
  {
    /// Has progress display been built?
    static const bool   display_progress =
  #ifdef _PIRANHA_DISPLAY_PROGRESS
      true;
#else
    false;
#endif
    /// Has integration with GSL libraries been enabled?
    static const bool   gsl_integration =
  #ifdef _PIRANHA_GSL
      true;
#else
    false;
#endif
    /// Is Piranha using Intel TBB for parallelism?
    static const bool   use_tbb =
  #ifdef _PIRANHA_TBB
      true;
#else
    false;
#endif
  };
}
#endif
