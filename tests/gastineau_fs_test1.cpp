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

using namespace piranha;
using namespace piranha::manipulators;

typedef dfs stype;

// Gastineau's test on Fourier series: power two of elp3 series without truncation, done 20 times.

int main()
{
	settings::set_debug(true);
	stype elp3("elp3.dfs");
	std::cout << "calculate pow 3"<<std::endl << std::flush;
	elp3 = elp3.pow(3);
	std::cout << "construct by copy"<<std::endl << std::flush;
	stype elp3a(elp3);
	std::cout << "square"<<std::endl << std::flush;
	truncators::Norm::set(1.0E10);
	elp3 * elp3a;
	std::cout << elp3.length() << '\n';
	std::cout << elp3a.length() << '\n';
}
