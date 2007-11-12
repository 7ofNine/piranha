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

#ifndef PIRANHA_COMMON_TYPEDEFS_H
#define PIRANHA_COMMON_TYPEDEFS_H

#include <boost/cstdint.hpp>
#include <complex>
#include <deque>
#include <string>
#include <valarray>
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
/// Maximum fast integer.
  typedef boost::int64_t max_fast_int;
#else
  typedef boost::int32_t max_fast_int;
#endif
/// Alias for the size of trigonometric containers.
  typedef uint16 trig_size_t;
/// Alias for complex doubles.
  typedef std::complex<double> complex_double;
/// Alias for array of int16.
  typedef std::valarray<int16> array_int16;
/// Alias for vector of int16.
  typedef std::vector<int16> vector_int16;
/// Alias for vector of double.
  typedef std::vector<double> vector_double;
/// Deque of strings.
  typedef std::deque<std::string> deque_string;
/// Alias for degree type.
  typedef int degree_type;
}

#endif
