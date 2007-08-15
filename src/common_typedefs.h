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

#include <complex>
#include <deque>
#include <string>
#include <valarray>
#include <vector>

// These are commonly used typedefs.

/// Alias for the size of trigonometric containers.
typedef unsigned short int trig_size_t;
/// Alias for the multipliers (elements of j vectors).
typedef short int mult_t;
/// Alias for complex doubles.
typedef std::complex<double> complex_double;
/// Alias for array of mult_t.
typedef std::valarray<mult_t> array_mult_t;
/// Alias for vector of mult_t.
typedef std::vector<mult_t> vector_mult_t;
/// Alias for vector of double.
typedef std::vector<double> vector_double;
/// Deque of strings. Convenience typedef.
typedef std::deque<std::string> deque_string;
/// Alias for exponent type.
typedef short int expo_type;
/// Alias for degree type.
typedef int degree_type;

#endif
