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

#include "../src/manipulators/zpoly.h"

// Exercise substitution.

using namespace piranha::manipulators;
using namespace piranha;

int main()
{
	psym x("x"), y("y"), z("z");
	zpoly f = zpoly(x) + zpoly(y) + zpoly(z), g = f.pow(max_fast_int(40));
	return !(g.sub(x,zpoly(x)+max_fast_int(1)).sub(x,zpoly(x)-max_fast_int(1)) == g);
}
