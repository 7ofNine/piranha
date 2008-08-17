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
settings::debug(true);
	dpoly x(psym("x")), y(psym("y")), z(psym("z")), t(psym("t")), u(psym("u"));
	max_fast_int one(1), two(2), three(3), five(5), twelve(12);
	dpoly f = (x + y + z*z*two + t*t*t*three + u.pow(five)*five + one).pow(twelve);
	dpoly g = (u + t + z*z*two + y*y*y*three + x.pow(five)*five + one).pow(twelve);
	f *= g;
	std::cout << f.length() << '\n';
	if (f.length() != 5821335) {
		return 1;
	}
	return 0;
}
