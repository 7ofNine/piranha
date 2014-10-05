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

#include "../src/manipulators/qps.h"

// Exercise special functions.

using namespace piranha::manipulators;
typedef qps ps;
using namespace piranha;

int main()
{
	int retval = 0;
	Psym x_("x");
	const ps x(x_);
	truncators::degree::set(80);
	retval += (x.besselJ(3) != -x.besselJ(-3));
	retval += (x.root(1) != x);
	retval += (x.pow(-3).root(-3) != x);
	retval += (x.sin() * x.sin() + x.cos() * x.cos() != 1);
	retval += (x.besselJ_div_m(-3,-3) != x.besselJ(-3) * x.pow(3));
	retval += (ps().besselJ(0) != 1 || ps().besselJ(1) != 0);
	// Test Bessel's differential equation.
	const int order = -5;
	const ps y = x.besselJ(order);
	retval += (
		x.pow(2)*y.partial("x",2)+x*y.partial("x")+
		(x.pow(2)-ps(order*order))*y != 0
	);
	// Relation between Jn and its derivatives.
	retval += (x.dbesselJ(5) * 2 != x.besselJ(4) - x.besselJ(6));
	return retval;
}
