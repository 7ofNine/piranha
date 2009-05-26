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
	// Expansion to order 400 of r/a in terms of e and M.
	int retval = 0;
	psym e("e"), M("M");
	degree_truncator::set(400);
	ps res(ps::r_a(ps(e),ps(M)));
	std::cout << res.length() << '\n';
	std::cout << res.atoms() << '\n';
	retval += (res.length() != 401 || res.atoms() != 80805);

	// Identity.
	degree_truncator::set(20);
	ps e_s(e), M_s(M);
	retval += (qps::EE(e_s,M_s).cos().sub(e,e_s-1).sub(e,e_s+1) != qps::cos_E(e_s,M_s));

	// Another identity.
	retval += (qps::r_a(e_s,M_s) * qps::sin_f(e_s,M_s) !=
		(-e_s*e_s + 1).root(2) * qps::sin_E(e_s,M_s));

	// Testing pow.
	retval += (ps::r_a(e_s,M_s).pow(3).root(3) != ps::r_a(e_s,M_s));

	// Trigonometric identity.
	retval += ((ps::cos_E(e_s,M_s).pow(2) + ps::sin_E(e_s,M_s).pow(2))
		!= 1);

	// r_a ** -2 == a_r ** 2 and vice-versa.
	retval += (ps::a_r(e_s,M_s).pow(2) != ps::r_a(e_s,M_s).pow(-2));
	retval += (ps::r_a(e_s,M_s).pow(2) != ps::a_r(e_s,M_s).pow(-2));

	// eipE tests.
	retval += (ps::eipE(e_s,M_s,5).abs() != 1);
	retval += (ps::eipE(e_s,M_s,5) * ps::eipE(e_s,M_s,-5) != 1);
	retval += (ps::eipE(e_s,M_s,1).pow(6) != ps::eipE(e_s,M_s,6));
	retval += (ps::eipE(e_s,M_s,1).pow(-1).pow(6) != ps::eipE(e_s,M_s,-6));

	return retval;
}
