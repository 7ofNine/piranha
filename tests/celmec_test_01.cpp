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

// Exercise celmec expansions, trigonometric functions, substitution and pow.

using namespace piranha::manipulators;
typedef qps ps;
using namespace piranha;

int main()
{
    settings::setMultiplicationAlgorithm(settings::ALGORITHM_PLAIN);
	// Expansion to order 401 of r/a in terms of e and M.
	int retval = 0;
    // todo: remove just fro some testing
    {
        Psym x("x"), y("y"), z("z");
        qps aSeries = 2*qps(x) + qps(y);
        qps bSeries =   qps(y) + qps(z);
        qps testSeries = aSeries + bSeries;

        Psym p("p"), q("q");
        qps testTrigSeries = 2*qps(p) + 3*qps(q);
        qps testResult = testSeries*(testTrigSeries.cos() + testTrigSeries.sin());
        
    }
    // todo end


	Psym e("e"), M("M");
	truncators::degree::set(4);
    std::vector< ps > flattened;
    
	    ps res(ps::r_a(ps(e), ps(M)));
	    std::cout << res.length() << '\n';
	    std::cout << res.atoms() << '\n';
	    retval += (res.length() != 401 || res.atoms() != 80805);
        flattened = res.flatten();
    
    std::cout<< "test 1: " << retval << '\n';
   
    Psym ep("e'"),Mp("M'");
    ps resp(ps::a_r(ps(ep), ps(Mp)));
    ps retp = res*resp;
//    retp.saveTo("testxxxx.qps");
 //   retp.print(std::cout);
//    return 0;
	// Identity.
	truncators::degree::set(20);
	ps e_s(e), M_s(M);
	retval += (qps::EE(e_s,M_s).cos().sub("e",e_s-1).sub("e",e_s+1) != qps::cos_E(e_s,M_s));
    std::cout<< "test 2: " << retval << '\n';

	// Another identity.
	retval += (qps::r_a(e_s,M_s) * qps::sin_f(e_s,M_s) !=
		(-e_s*e_s + 1).root(2) * qps::sin_E(e_s,M_s));
    std::cout<< "test 3: " << retval << '\n';
	// Testing pow.
	retval += (ps::r_a(e_s,M_s).pow(3).root(3) != ps::r_a(e_s,M_s));
    std::cout<< "test 4: " << retval << '\n';
	// Trigonometric identity.
	retval += ((ps::cos_E(e_s,M_s).pow(2) + ps::sin_E(e_s,M_s).pow(2))
		!= 1);
    std::cout<< "test 5: " << retval << '\n';
	// r_a ** -2 == a_r ** 2 and vice-versa.
	retval += (ps::a_r(e_s,M_s).pow(2) != ps::r_a(e_s,M_s).pow(-2));
     std::cout<< "test 6: " << retval << '\n';
	retval += (ps::r_a(e_s,M_s).pow(2) != ps::a_r(e_s,M_s).pow(-2));
     std::cout<< "test 7: " << retval << '\n';

	// eipE tests.
	retval += (ps::eipE(e_s,M_s,5).abs() != 1);
     std::cout<< "test 8: " << retval << '\n';
	retval += (ps::eipE(e_s,M_s,5) * ps::eipE(e_s,M_s,-5) != 1);
     std::cout<< "test 9: " << retval << '\n';
	retval += (ps::eipE(e_s,M_s,1).pow(6) != ps::eipE(e_s,M_s,6));
     std::cout<< "test 10: " << retval << '\n';
	retval += (ps::eipE(e_s,M_s,1).pow(-1).pow(6) != ps::eipE(e_s,M_s,-6));
     std::cout<< "test 11: " << retval << '\n';
    std::cout << retval << std::endl << std::flush;
	return retval;
}
