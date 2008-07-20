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

#include <iostream>

#include "../src/manipulators/dfs.h"

using namespace piranha;
using namespace piranha::manipulators;

typedef dfs stype;

// Gastineau's test on Fourier series: power two of elp3 series without truncation, done 20 times.

int main()
{
  base_norm_truncator::set(0);
  stype elp3("elp3.fs"), elp3a(elp3);
  for (size_t i = 0; i < 20; ++i)
  {
    std::cout << (elp3*elp3a).length() << std::endl;
  }
  return 0;
}
