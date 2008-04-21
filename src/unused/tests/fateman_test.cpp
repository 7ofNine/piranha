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

#include "../src/piranha.h"

typedef piranha::ipoly<long long int,size_t,unsigned short int> poly_type;
typedef poly_type::const_iterator const_iterator;

// Compute q*(q+1), where q=(1+x+y+z)^20. Do it 10 times.
int main()
{
  for (size_t i=0;i<10;++i)
  {
    poly_type::vector_expo v(3);
    poly_type p1(1,v);
    v[0]=1;
    poly_type p2(1,v);
    v[0]=0;
    v[1]=1;
    poly_type p3(1,v);
    v[1]=0;
    v[2]=1;
    poly_type p4(1,v);
    p1+=p2+=p3+=p4;
    poly_type q(p1);
    std::cout << "Here comes p1:\n" << p1;
    //std::cout << std::scientific << std::setprecision(15) << "p1 evals to: " << eval_poly(p1) << '\n';
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    p1*=q;
    //std::cout << std::scientific << std::setprecision(15) << "p1 evals to: " << eval_poly(p1) << '\n';
    poly_type r(p1);
    r+=poly_type(1,poly_type::vector_expo(3));
    p1*=r;
    //std::cout << std::scientific << std::setprecision(15) << "p1 evals to: " << eval_poly(p1) << '\n';
    std::cout << "Degree is:" << p1.g_degree() << '\n';
  }
  return 0;
}
