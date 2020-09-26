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

#include "../src/manipulators/dpoly.h"

// Small polynomials multiplication test.
// s*(s+1)
// where s = (x+1)^20. Do it 10^4 times.

using namespace piranha;
using namespace piranha::manipulators;
typedef dpoly poly;

int main()
{
	poly a_("small_poly_test.dpoly"), b_(a_);
	for (size_t i = 0; i < 10000; ++i)
	{
		poly a(a_), b(b_);
		a = a.pow((MaxFastInt)20);
		poly c(a);
		c+=(MaxFastInt)1;
		a*=c;
	}
	return 0;
}
