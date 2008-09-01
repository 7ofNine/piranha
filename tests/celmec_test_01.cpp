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
	base_degree_truncator::set(400);
	ps res(ps::r_a(ps(e),ps(M)));
	std::cout << res.length() << '\n';
	std::cout << res.atoms() << '\n';
	retval += (res.length() != 401 || res.atoms() != 80805);

	// Identity.
	base_degree_truncator::set(20);
	ps e_s(e), M_s(M);
	retval += (qps::EE(e,M).cos().sub(e,e_s-max_fast_int(1)).sub(e,e_s+max_fast_int(1)) != qps::cos_E(e,M));

	// Testing pow.
	retval += (ps::r_a(ps(e),ps(M)).pow(max_fast_int(3)).root(max_fast_int(3)) != ps::r_a(ps(e),ps(M)));

	return retval;
}
