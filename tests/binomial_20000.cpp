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

#include "../src/manipulators/dpoly.h"

// Calculate:
// ((x+y)**2)*10000

using namespace piranha::manipulators;
typedef dpoly poly;
using namespace piranha;

int main()
{
settings::set_debug(false);
//settings::set_nthread(16); // temporary: remove
try{
	poly x(Psym("x")), y(Psym("y"));
	poly res((x+y).pow(2).pow(10000));
	if (res.length() != 20001) {
		return 1;
	}
	std::cout << res.length() << '\n';
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
	return 0;
}
