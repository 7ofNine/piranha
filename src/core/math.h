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

#ifndef PIRANHA_MATH_H
#define PIRANHA_MATH_H

#include <boost/static_assert.hpp>

namespace piranha
{
  /// Meta-programmed functor for the calculation of base-2 logarithm.
  /**
   * Result is retrieved through the lg::value const member function.
   */
  template <int N>
    struct lg
  {
    BOOST_STATIC_ASSERT(N > 0 and (N % 2) == 0);
    /// Value of the base-2 logarithm of N.
    static const size_t value = lg<(N >> 1)>::value+1;
  };

  template <>
    struct lg<1>
  {
    static const size_t value = 0;
  };
}

#endif
