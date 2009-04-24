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
#include "../src/manipulators/zpoly.h"

// Exercise substitution.

using namespace piranha::manipulators;
using namespace piranha;

int main()
{
	{
		// This one fails if we do not handle correctly args_tuple inside substitution.
		qps e(psym("e")), ph(psym("ph")), th(psym("th"));
		degree_truncator::set(10);
		(e*th.cos()+1).pow(-1).sub(psym::get("th"),ph.pow(2));
	}

	int retval = 0;
	psym x("x"), y("y"), z("z");
	zpoly f = zpoly(x) + zpoly(y) + zpoly(z), g = f.pow(40);
	retval += !(g.sub(x,zpoly(x)+1).sub(x,zpoly(x)-1) == g);
	retval += (qps(x).cos().ei_sub(x,qpsc(std::complex<double>(1,0))) != 1);
	retval += (qps(x).sin().ei_sub(x,qpsc(std::complex<double>(1,0))) != 0);
	return retval;
}
