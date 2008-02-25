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

#ifndef PIRANHA_COMMON_TYPEDEFS_H
#define PIRANHA_COMMON_TYPEDEFS_H

#include <boost/cstdint.hpp>
#include <complex>
#include <deque>
#include <string>
#include <vector>

namespace piranha
{
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
#ifdef _PIRANHA_64BIT
  /// Maximum fast integer (64-bits).
  typedef boost::int64_t max_fast_int;
  /// Maximum fast unsigned integer (64-bits).
  typedef boost::uint64_t max_fast_uint;
#else
  /// Maximum fast integer (32-bits).
  typedef boost::int32_t max_fast_int;
  /// Maximum fast unsigned integer (32-bits).
  typedef boost::uint32_t max_fast_uint;
#endif
  /// Alias for complex doubles.
  typedef std::complex<double> complex_double;
  /// Alias for vector of double.
  typedef std::vector<double> vector_double;
  /// Deque of strings.
  typedef std::deque<std::string> deque_string;
  /// Layout element, to be used in series merging.
  typedef std::pair<bool,size_t> layout_element;
  /// Layout type, to be used in series merging.
  typedef std::vector<layout_element> layout_type;
}
#endif
