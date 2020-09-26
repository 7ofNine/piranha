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

//#include "../src/manipulators/dpoly.h"
#include "piranha.h"
// Pearce's sparse polynomial multiplication test 2. Calculate:
// f^5*g^5
// where 
// f = x1*x2+x1*x10+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7+x7*x8+x8*x9+x9*x10+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+1;
// g = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2+x9^2+x10^2+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+1;

using namespace piranha::manipulators;
typedef dpoly poly;

int main()
{
  try{

    dpoly t1("pearce_sparse3.dpoly");
    dpoly t2("pearce_sparse4.dpoly");
    dpoly u1(t1);
    dpoly u2(t2);
    t1*=u1;
    t1*=u1;
    t1*=u1;
    t1*=u1;

    t2*=u2;
    t2*=u2;
    t2*=u2;
    t2*=u2;

    t1*=t2;

    std::cout << std::endl << t1.length() << std::endl;
 } catch(std::exception ex)
 {
    std::cout << ex.what() << std::endl;
    return 1;
 } 

  return 0;
}
