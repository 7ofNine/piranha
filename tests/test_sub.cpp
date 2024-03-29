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

//#include "../src/manipulators/qps.h"
//#include "../src/manipulators/zpoly.h"
#include "piranha.h"

// Exercise substitution.

using namespace piranha::manipulators;
using namespace piranha;

int main()
{
	{
		// This one fails if we do not handle correctly argsTuple inside substitution.
		qps e(Psym("e"));
		qps ph(Psym("ph"));
		qps th(Psym("th"));
		truncators::Degree::set(10); // this is global
		qps trigsub = (e*th.cos()+1).pow(-1).sub("th",ph.pow(2));
	}
    truncators::Degree::unset(); // deactivate truncation
	int retval = 0;
	Psym x("x");
	Psym y("y");
	Psym z("z");

    zpoly f = zpoly(x) + zpoly(y) + zpoly(z);
    zpoly g = f.pow(40); ///original 40

    zpoly tempg = g.sub("x", zpoly(x) + 1).sub("x", zpoly(x) - 1);
    tempg.printPlain(std::cout);
	retval += !(tempg == g);
	
    retval += (qps(x).cos().eiSubstitute("x", qpsc(std::complex<double>(1, 0))) != 1);
	retval += (qps(x).sin().eiSubstitute("x", qpsc(std::complex<double>(1, 0))) != 0);
    
    std::cout << retval << std::endl;
	return retval;
}
