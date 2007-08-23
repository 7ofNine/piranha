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

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include "base_pseries.h"
#include "double_cf.h"
#include "polynomial.h"
#include "simple_term.h"
#include "trig_array.h"

// Toolboxes.
#include "operators_toolbox.h"
#include "symbol_limiting_elementary_math_toolbox.h"
#include "differential_toolbox.h"

namespace piranha
{
/// Extract minimum exponent from terms with polynomials as coefficients.
  template <class Term>
    struct min_expo_extractor
  {
    typedef expo_type result_type;
    expo_type operator()(const Term &t) const
    {
      return t.g_cf()->g_min_expo();
    }
  };

/// Extract degree from terms with polynomials as coefficients.
  template <class Term>
    struct degree_extractor
  {
    typedef degree_type result_type;
    degree_type operator()(const Term &t) const
    {
      return t.g_cf()->g_degree();
    }
  };

/// Degree-based indices for base_pseries.
/**
 * This class specifies the following indices to be used in piranha::base_pseries: a hashed index for the
 * identification
 * of terms and a minimum exponent sorted index to discard terms in multiplications. The class is to be used as the I
 * parameter in piranha::base_pseries classes.
 */
  template <class Cf, class Trig, template <class, class> class Term>
    struct degree_based_index
  {
    typedef boost::multi_index::indexed_by <
      boost::multi_index::ordered_unique <
      boost::multi_index::composite_key <
      Term<Cf, Trig>,
      min_expo_extractor<Term<Cf, Trig> >,
      degree_extractor<Term<Cf, Trig> >,
      boost::multi_index::const_mem_fun < Term<Cf, Trig>, const bool &,
      &Term<Cf, Trig>::g_flavour > ,
      boost::multi_index::const_mem_fun < Term<Cf, Trig>, const Trig &,
      &Term<Cf, Trig>::g_trig_ref >
      >
      >,
      boost::multi_index::hashed_unique <
      boost::multi_index::composite_key <
      Term<Cf, Trig>,
      boost::multi_index::const_mem_fun < Term<Cf, Trig>, const bool &,
      &Term<Cf, Trig>::g_flavour > ,
      boost::multi_index::const_mem_fun < Term<Cf, Trig>, const Trig &,
      &Term<Cf, Trig>::g_trig_ref >
      >
      >
      > type;
  };

/// Derived class for symbolic Poisson series.
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I>
    class sps:
  public base_pseries<Cf,Trig,Term,I,sps<Cf,Trig,Term,I> >,
    public symbol_limiting_elementary_math_toolbox<sps<Cf,Trig,Term,I> >,
    public real_operators_toolbox<sps<Cf,Trig,Term,I> >,
    public differential_toolbox<sps<Cf,Trig,Term,I> >
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
      explicit sps(const cf_type &c, const sps &model):ancestor::base_pseries(c,model)
        {}
/// Constructor from int.
      explicit sps(int n)
      {
        ancestor::generic_builder(n);
      }
/// Constructor from double.
      explicit sps(const double &x)
      {
        ancestor::generic_builder(x);
      }
/// Constructor from psymbol.
      explicit sps(const psymbol &psym, psymbol::type ptype):ancestor::base_pseries(psym,ptype)
        {}
  };

  typedef sps<polynomial<double_cf>,trig_array,simple_term,degree_based_index> sp;
}
#endif
