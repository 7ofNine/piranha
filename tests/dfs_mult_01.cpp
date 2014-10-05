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

#include <iostream>

#include "../src/manipulators/dfs.h"
#include "../src/manipulators/qps.h"

using namespace piranha;
using namespace piranha::manipulators;

typedef dfs stype;


// Double-precision fs multiplication test.

int main()
{
	settings::set_debug(true);
//	stype elp3("elp3.dfs");
//	elp3 *= elp3;
//	std::cout << (elp3*elp3).length() << std::endl;
    qps test("testxxxx.qps");
    test *= test;
	return 0;
}
