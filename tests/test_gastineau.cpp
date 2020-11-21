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

//#include "../src/manipulators/dfs.h"
#include "piranha.h"

using namespace piranha;
using namespace piranha::manipulators;

typedef dfs stype;

// Gastineau's test on Fourier series: power two of elp3 series without truncation, done 20 times.

int main()
{
    int retval = 0;
#ifdef DEBUG
	settings::set_debug(false);
#endif
try{
	stype elp3("elp3.dfs");
	std::cout << "input :" << " length: " << elp3.length() << " atoms: " << elp3.atoms() << std::endl;
    retval +=(elp3.length() != 702)||(elp3.atoms() != 1404);
    std::cout << "read: " << retval;

    elp3 = elp3.pow(3);
	stype elp3a(elp3);
    std::cout << "copy construct :" << " length: " << elp3a.length() << " atoms: " << elp3a.atoms() << std::endl;
    retval += (elp3a.length() != elp3.length())||(elp3a.length() != 60205);   // currently elp3.length result in 60204 ???
    std::cout << "copy construct : " << retval << std::endl;

	dfs square  = elp3 * elp3a;
    retval += ((square.length() != 980359)) || ((square.atoms() != 1960718));
   	std::cout << "square : " << " length :" << square.length() << " atoms: " << square.atoms() << std::endl;
	std::cout << "gastineau_fs_test1: " << retval<< std::endl;
}catch(p_base_exception &pex)
    {
	    std::cout << pex.what() << std::endl;
	    return -1;
    }
catch(std::exception &ex)
    {
        std::cout << ex.what() << std::endl;
        return -1;
    }

    return retval;
}
