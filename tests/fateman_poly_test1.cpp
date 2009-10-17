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

#include "../src/manipulators/qqpoly.h"

// Fateman's polynomial multiplication test number 1. Calculate:
// s*(s+1)
// where s = (1+x+y+z+t)^30

using namespace piranha;
using namespace piranha::manipulators;
typedef qqpoly poly;

int main()
{
//settings::debug(true);
//settings::set_nthread(1);
  psym x("x"), y("y"), z("z"), t("t");
  poly a = poly(x) + poly(y) + poly(z) + poly(t) + poly(1) , b(a);
  a = std::pow(a,15);
  poly c(a);
  c+=1;
  a*=c;

  std::cout << std::endl << a.length() << std::endl;

  return 0;
}
