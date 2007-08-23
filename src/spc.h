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

#ifndef PIRANHA_SPC_H
#define PIRANHA_SPC_H

#include "complex_toolbox.h"
#include "sp.h"

namespace std
{
/// Complex specialization for sp class.
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I>
    struct complex<piranha::sps<Cf,Trig,Term,I> >:
    public piranha::base_pseries <complex<Cf>,Trig,Term,I,complex<piranha::sps<Cf,Trig,Term,I> > >,
    public piranha::symbol_limiting_elementary_math_toolbox<sps<Cf,Trig,Term,I> >,
    public piranha::complex_operators_toolbox<piranha::sps<Cf,Trig,Term,I> >,
    public piranha::complex_toolbox<piranha::sps<Cf,Trig,Term,I> >
  {
    public:
/// Alias for the ancestor.
      typedef piranha::base_pseries<complex<Cf>,Trig,Term,I,complex<piranha::sps<Cf,Trig,Term,I> > > ancestor;
/// Alias for term type.
      typedef typename ancestor::term_type term_type;
/// Alias for the iterator to sorted index.
      typedef typename ancestor::it_s_index it_s_index;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Alias for real counterpart.
      typedef piranha::sps<Cf,Trig,Term,I> real_ps;
/// Alias for real coefficient.
      typedef Cf real_cf_type;
/// Alias for self.
      typedef complex<piranha::sps<Cf,Trig,Term,I> > self;
/// Default constructor.
      complex():ancestor::base_pseries()
        {}
/// Copy constructor.
      complex(const complex &p):ancestor::base_pseries(p)
        {}
/// Constructor from filename.
      explicit complex(const std::string &filename):ancestor::base_pseries(filename)
        {}
/// Constructor from coefficient.
      explicit complex(const cf_type &c, const complex &model):ancestor::base_pseries(c,model)
        {}
/// Constructor from pair of real coefficients.
/*explicit complex(const real_cf_type &a, const real_cf_type &b):
  ancestor::base_pseries(cf_type(a,b))
  {}*/
/// Constructor from complex.
      explicit complex(const complex_double &c)
      {
        ancestor::generic_builder(c);
      }
/// Constructor from pair of doubles.
      explicit complex(const double &a, const double &b)
      {
        ancestor::generic_builder(complex_double(a,b));
      }
/// Constructor from real series.
// FIXME: here and below we are discarding lin_args.
// TODO: can we re-use some function from complex_toolbox to achieve this result?
// If so, ditch the term_type typedef which is used just here. Also the iterator typedef..
      explicit complex(const real_ps &p)
      {
        p_assert(ancestor::merge_args(p));
        term_type term;
        typename real_ps::it_s_index it=p.s_index().begin(), it_f=p.s_index().end();
        for (;it!=it_f;++it)
        {
          *term.s_cf()=cf_type(*it->g_cf());
          *term.s_trig()=*it->g_trig();
          term.s_flavour()=it->g_flavour();
          ancestor::insert(term);
        }
      }
/// Constructor from real and imaginary series.
      explicit complex(const real_ps &p, const real_ps &q)
      {
        build_from_components(p,q);
      }
/// Constructor from real and imaginary series from filenames.
      explicit complex(const std::string &file1, const std::string &file2)
      {
        build_from_components(real_ps(file1),real_ps(file2));
      }
  };
}

namespace piranha
{
  typedef std::complex<sp> spc;
}
#endif
