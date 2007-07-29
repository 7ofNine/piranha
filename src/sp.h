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

#ifndef PIRANHA_SP_H
#define PIRANHA_SP_H

#include "base_pseries.h"
#include "double_cf.h"
#include "polynomial.h"
#include "ps_term.h"
#include "trig_array.h"

namespace piranha
{
/// Derived class for symbolic Poisson series.
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I>
    class sps:
    public base_pseries<Cf,Trig,Term,I,sps<Cf,Trig,Term,I> >
  {
    public:
/// Alias for parent class.
      typedef piranha::base_pseries<Cf, Trig, Term, I, sps<Cf,Trig,Term,I> > ancestor;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Alias for term type.
      typedef typename ancestor::term_type term_type;
/// Alias for self.
      typedef sps<Cf,Trig,Term,I> self;
// Ctors.
      sps():ancestor::base_pseries()
        {}
      sps(const sps &p):ancestor::base_pseries(p)
        {}
      explicit sps(const std::string &filename):ancestor::base_pseries(filename)
        {}
      explicit sps(const cf_type &c):ancestor::base_pseries(c)
        {}
/// Constructor from int.
      explicit sps(int n):ancestor::base_pseries(cf_type(n))
        {}
/// Constructor from double.
      explicit sps(const double &x):ancestor::base_pseries(cf_type(x))
        {}
/// Constructor from psymbol.
      explicit sps(const psymbol &psym, psymbol::type ptype):ancestor::base_pseries(psym,ptype)
        {}
  }
  ;


  typedef sps<polynomial<double_cf>,trig_array,ps_term,norm_based_index> sp;
}


#endif
