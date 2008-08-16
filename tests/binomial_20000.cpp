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

#include "../src/manipulators/zpoly.h"

// Calculate:
// ((x+y)**2)*10000

using namespace piranha::manipulators;
typedef zpoly poly;
using namespace piranha;

int main()
{
	settings::debug(true);
	settings::memory_limit(1000000000);
	poly x(psym("x")), y(psym("y"));
	poly res((x+y).pow(max_fast_int(2)).pow(max_fast_int(10000)));
	if (res.length() != 20001) {
		return 1;
	}
	std::cout << res.length() << '\n';
	return 0;
}
