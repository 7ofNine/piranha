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

#include "bits/packed_int_array.h"

namespace piranha
{
  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size64;

  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size32;

  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size16;

  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::size8;

#ifdef _PIRANHA_SSE2
  template <int Dim, int Bits>
    const uint8 packed_int_array<Dim,Bits>::m128n;
#endif
}
