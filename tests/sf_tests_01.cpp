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
	psym x_("x");
	const ps x(x_);
	base_degree_truncator::set(80);
	retval += (x.besselJ(max_fast_int(3)) != -x.besselJ(max_fast_int(-3)));
	retval += (x.root(max_fast_int(1)) != x);
	retval += (x.pow(max_fast_int(-3)).root(max_fast_int(-3)) != x);
	retval += (x.sin() * x.sin() + x.cos() * x.cos() != max_fast_int(1));
	retval += (x.pow(max_fast_int(10)) * x.besselJ_div_m(max_fast_int(-6),max_fast_int(10)) != x.besselJ(max_fast_int(-6)));
	retval += (x.besselJ_div_m(max_fast_int(-3),max_fast_int(-3)) != x.besselJ(max_fast_int(-3)) * x.pow(max_fast_int(3)));
	retval += (ps().besselJ(max_fast_int(0)) != max_fast_int(1) || ps().besselJ(max_fast_int(1)) != max_fast_int(0));
	// Test Bessel's differential equation.
	const max_fast_int order = -5;
	const ps y = x.besselJ(order);
	retval += (
		x.pow(max_fast_int(2))*y.partial(x_,max_fast_int(2))+x*y.partial(x_)+
		(x.pow(max_fast_int(2))-ps(order*order))*y != max_fast_int(0)
	);
	// Relation between Jn and its derivatives.
	retval += (x.dbesselJ(5) * max_fast_int(2) != x.besselJ(4) - x.besselJ(6));
	return retval;
}
