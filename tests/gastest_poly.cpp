/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

typedef piranha::ipoly<double,size_t,unsigned short int> poly_type;
typedef poly_type::const_iterator const_iterator;

// Evaluate polynomial with fixed value 1 for all variables.
inline double eval_poly(const poly_type &p)
{
  poly_type::vector_expo ev((size_t)(p.g_width()));
  std::vector<double> values((size_t)(p.g_width()));
  for (size_t i=0;i<p.g_width();++i)
  {
    values[i]=1;
  }
  double retval=0, tmp;
  for (const_iterator it=p.begin();it!=p.end();++it)
  {
    p.decode(it->g_index(),ev);
    tmp=it->g_cf();
    for (size_t i=0;i<p.g_width();++i)
    {
      tmp*=std::pow(values[i],ev[i]);
    }
    retval+=tmp;
  }
  return retval;
}

// Compute s*(s+1), where s=(1+x+y+t+u)^14. Do it 10 times.
int main(int argc, char **argv)
{
  for (size_t i=0;i<10;++i)
  {
    poly_type::vector_expo v(4);
    poly_type p1(1.,v);
    v[0]=1;
    poly_type p2(1.,v);
    v[0]=0;
    v[1]=1;
    poly_type p3(1.,v);
    v[1]=0;
    v[2]=1;
    poly_type p4(1.,v);
    v[2]=0;
    v[3]=1;
    poly_type p5(1.,v);
    p1+=p2+=p3+=p4+=p5;
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
    //std::cout << std::scientific << std::setprecision(15) << "p1 evals to: " << eval_poly(p1) << '\n';
    poly_type r(p1);
    r+=poly_type(1.,poly_type::vector_expo(4));
    p1*=r;
    //std::cout << std::scientific << std::setprecision(15) << "p1 evals to: " << eval_poly(p1) << '\n';
    std::cout << "Degree is:" << p1.g_degree() << '\n';
  }
  return 0;
}
