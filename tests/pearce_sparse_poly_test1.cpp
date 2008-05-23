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

#include "../src/piranha.h"

// Pearce's sparse polynomial multiplication test 1. Calculate:
// f*g
// where 
// f = (1 + x + y + 2*z^2 + 3*t^3 + 5*u^5)^12:
// g = (1 + u + t + 2*z^2 + 3*y^3 + 5*x^5)^12:

using namespace piranha::manipulators;
typedef dpoly poly;
using namespace piranha;

int main()
{
  dpoly t1("pearce_sparse1.dpoly");
  dpoly t2("pearce_sparse2.dpoly");
  dpoly u1(t1);
  dpoly u2(t2);
  t1 = std::pow(t1,(max_fast_int)12);
  t2 = std::pow(t2,(max_fast_int)12);
  t1 *= t2;

  std::cout << std::endl << t1.nth_index<0>().size() << std::endl;

  return 0;
}
