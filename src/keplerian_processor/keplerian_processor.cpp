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

#include "keplerian_processor.h"

namespace piranha
{
  gsp keplerian_processor::E_pow_e(const psymbol &e, const psymbol &M, int n)
    {
      if (n<0)
      {
        std::cout << "WARNING: a positive integer must be supplied for the order of the development." << std::endl;
        std::cout << "WARNING: returning empty series." << std::endl;
        return gsp();
      }
      gsp retval(M,psymbol::trig);
      for (int i=1;i<=n;++i)
      {
        gsp tmp=math::pow_besselJ<gsp,mpz_class>(i,gsp(e,psymbol::cf)*i,(unsigned int)n);
        //retval+=(math::pow_besselJ(i,Gsp.gsp(e,psymbol::cf)*i,6)*2)/n*(Gsp.gsp(psymbol("M",0,1),trig)*n).sine()
      }
      return retval;
    }
}
