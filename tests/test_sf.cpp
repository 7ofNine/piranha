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
#include "piranha.h"

// Exercise special functions.

using namespace piranha::manipulators;
typedef qps ps;
using namespace piranha;

int main()
{
	int retval = 0;
	Psym x_("x");
	const ps x(x_);
	truncators::Degree::set(80);
    ps lhs;
    ps rhs;
    lhs = x.besselJ(3);
    rhs = -x.besselJ(-3);
	//retval += (x.besselJ(3) != -x.besselJ(-3));
    retval +=(lhs != rhs);
    std::cout << "test 1: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 1: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 1: " << retval << std::endl;

    lhs = x.root(1);
    rhs = x;
    //retval += (x.root(1) != x);
    retval +=(lhs != rhs);
    std::cout << "test 2: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 2: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 2: " << retval << std::endl;

    lhs = x.pow(-3).root(-3);
    rhs = x;
    //retval += (x.pow(-3).root(-3) != x);
    retval +=(lhs != rhs);
    std::cout << "test 3: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 3: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 3: " << retval << std::endl;

    lhs = x.sin() * x.sin(); 
    rhs =  x.cos() * x.cos();
	//retval += (x.sin() * x.sin() + x.cos() * x.cos() != 1);
    retval += (lhs + rhs != 1);
    std::cout << "test 4: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 4: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 4: " << retval << std::endl;

    lhs = x.besselJ_div_m(-3,-3);
    rhs = x.besselJ(-3) * x.pow(3);
	//retval += (x.besselJ_div_m(-3,-3) != x.besselJ(-3) * x.pow(3));
    retval +=(lhs != rhs); 
    std::cout << "test 5: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 5: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 5: " << retval << std::endl;


    lhs = ps().besselJ(0);
    rhs = ps().besselJ(1);
	//retval += (ps().besselJ(0) != 1 || ps().besselJ(1) != 0);
    retval +=((lhs != 1) || (rhs != 0)); 
    std::cout << "test 6: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 6: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 6: " << retval << std::endl;

	// Test Bessel's differential equation.
	const int order = -5;
	const ps y = x.besselJ(order);

    lhs = x.pow(2)*y.partial("x", 2)+x*y.partial("x");
    rhs = (x.pow(2)-ps(order*order))*y;
	//retval += (x.pow(2)*y.partial("x", 2)+x*y.partial("x") + (x.pow(2)-ps(order*order))*y != 0);
    retval += ((lhs + rhs) != 0);
    std::cout << "test 7: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 7: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 7: " << retval << std::endl;

    lhs = x.dbesselJ(5) * 2;
    rhs = x.besselJ(4) - x.besselJ(6);
	// Relation between Jn and its derivatives.
	//retval += (x.dbesselJ(5) * 2 != x.besselJ(4) - x.besselJ(6));
    retval += (lhs != rhs);
    std::cout << "test 8: " << "lhs.length: " << lhs.length() << " lhs.atoms(): " << lhs.atoms() << std::endl;
    std::cout << "test 8: " << "rhs.length: " << rhs.length() << " rhs.atoms(): " << rhs.atoms() << std::endl;
    std::cout << "test 8: " << retval << std::endl;
	return retval;
}
