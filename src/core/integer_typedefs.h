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

#ifndef PIRANHA_INTEGER_TYPEDEFS_H
#define PIRANHA_INTEGER_TYPEDEFS_H

#include <boost/cstdint.hpp>

namespace piranha
{
  /// Integer selector.
  /**
   * It detects wehther the platform is 64bit or 32bit, and sets maximum and minimum
   * "fast" integer types accordingly, using the boost integer libraries. If the platform is
   * other than 32bit or 64bit it won't define any type.
   */
  template <int SizeOfPointer>
    struct int_selector {};

  // Specialization for 32bit archs.
  template <>
    struct int_selector<4>
  {
    typedef boost::int32_t  max_fast_int;
    typedef boost::uint32_t max_fast_uint;
  };

  // Specialization for 64bit archs.
  template <>
    struct int_selector<8>
  {
    typedef boost::int64_t  max_fast_int;
    typedef boost::uint64_t max_fast_uint;
  };

  // These are commonly used typedefs.
  /// Alias for 8bit integer.
  typedef boost::int8_t int8;
  /// Alias for 8bit unsigned integer.
  typedef boost::uint8_t uint8;
  /// Alias for 16bit integer.
  typedef boost::int16_t int16;
  /// Alias for unsigned 16bit integer.
  typedef boost::uint16_t uint16;
  /// Alias for 32bit integer.
  typedef boost::int32_t int32;
  /// Alias for 64bit integer.
  typedef boost::int64_t int64;
  /// Maximum fast integer, detected through piranha::int_selector.
  typedef int_selector<sizeof(void *)>::max_fast_int max_fast_int;
  /// Maximum fast unsigned integer, detected through piranha::int_selector.
  typedef int_selector<sizeof(void *)>::max_fast_uint max_fast_uint;
}
#endif
