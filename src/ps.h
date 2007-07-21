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

#ifndef PIRANHA_PS_H
#define PIRANHA_PS_H

#include <complex>

#include "base_pseries.h"
#include "operators_toolbox.h"
#include "norm_based_elementary_math_toolbox.h"
#include "math_toolbox.h"
#include "complex_toolbox.h"
#include "differential_toolbox.h"
#include "trigonometric_toolbox.h"

namespace piranha
{
/// Default derived class for Poisson series.
/**
 * This derived class provides support for complex arithmetic and for trigonometric operations on
 * series.
 */
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I>
    class ps:
  public base_pseries<Cf,Trig,Term,I,ps<Cf,Trig,Term,I> >,
    public real_operators_toolbox<ps<Cf,Trig,Term,I> >,
    public math_toolbox<ps<Cf,Trig,Term,I> >,
    public common_trig_toolbox<ps<Cf,Trig,Term,I>,ps<Cf,Trig,Term,I> >,
    public real_trig_toolbox<ps<Cf,Trig,Term,I> >,
    public differential_toolbox<ps<Cf,Trig,Term,I> >,
    public norm_based_elementary_math_toolbox<ps<Cf,Trig,Term,I> >
  {
    public:
      typedef piranha::base_pseries<Cf, Trig, Term, I, ps<Cf,Trig,Term,I> > ancestor;
      typedef typename ancestor::term_type term_type;
      typedef typename ancestor::it_s_index it_s_index;
      typedef typename ancestor::r_it_s_index r_it_s_index;
      typedef typename ancestor::cf_type cf_type;
/// Alias for self.
      typedef ps real_ps;
/// Alias for complex coefficient type.
      typedef std::complex<cf_type> complex_cf_type;
/// Alias for complex counterpart.
      typedef std::complex<ps<cf_type,Trig,Term,I> > complex_ps;
// Ctors.
      ps():ancestor::base_pseries()
        {}
      ps(const ps &p):ancestor::base_pseries(p)
        {}
      explicit ps(const std::string &filename):ancestor::base_pseries(filename)
        {}
      explicit ps(const cf_type &c):ancestor::base_pseries(c)
        {}
/// Constructor from int.
      explicit ps(int n):ancestor::base_pseries(cf_type(n))
        {}
/// Constructor from double.
      explicit ps(const double &x):ancestor::base_pseries(cf_type(x))
        {}
/// Constructor from psymbol.
      explicit ps(const psymbol &psym, psymbol::type ptype):ancestor::base_pseries(psym,ptype)
        {}
  }
  ;
}


namespace std
{
// COMPLEX COUNTERPART
/// Complex specialization for default derived class.
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I>
    struct complex<piranha::ps<Cf,Trig,Term,I> >:
  public piranha::base_pseries <complex<Cf>,Trig,Term,I,complex<piranha::ps<Cf,Trig,Term,I> > >,
    public piranha::complex_operators_toolbox<piranha::ps<Cf,Trig,Term,I> >,
    public piranha::common_trig_toolbox<complex<piranha::ps<Cf,Trig,Term,I> >,piranha::ps<Cf,Trig,Term,I> >,
    public piranha::complex_toolbox<piranha::ps<Cf,Trig,Term,I> >,
    public piranha::differential_toolbox<complex<piranha::ps<Cf,Trig,Term,I> > >,
    public piranha::norm_based_elementary_math_toolbox<complex<piranha::ps<Cf,Trig,Term,I> > >
  {
    public:
/// Alias for the ancestor.
      typedef piranha::base_pseries<complex<Cf>,Trig,Term,I,complex<piranha::ps<Cf,Trig,Term,I> > > ancestor;
/// Alias for term type.
      typedef typename ancestor::term_type term_type;
/// Alias for the iterator to sorted index.
      typedef typename ancestor::it_s_index it_s_index;
/// Alias for the iterator to the reversed sorted index.
      typedef typename ancestor::r_it_s_index r_it_s_index;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Alias for self.
      typedef complex complex_ps;
/// Alias for real counterpart.
      typedef piranha::ps<Cf,Trig,Term,I> real_ps;
/// Alias for real coefficient.
      typedef Cf real_cf_type;
// TODO: this is shared with above, find a way to share.
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
      explicit complex(const cf_type &c):ancestor::base_pseries(c)
        {}
/// Constructor from pair of real coefficients.
      explicit complex(const real_cf_type &a, const real_cf_type &b):ancestor::base_pseries(cf_type(a,b))
        {}
/// Constructor from complex.
      explicit complex(const complex_double &c):ancestor::base_pseries(cf_type(c))
        {}
/// Constructor from pair of doubles.
      explicit complex(const double &a, const double &b):ancestor::base_pseries(cf_type(complex_double(a,b)))
        {}
/// Constructor from real series.
// FIXME: here and below we are discarding lin_args.
// TODO: can we re-use some function from complex_toolbox to achieve this result?
      explicit complex(const real_ps &p):ancestor::base_pseries()
      {
        p_assert(ancestor::merge_args(p));
        term_type term;
        typename real_ps::it_s_index it=p.s_index().begin(), it_f=p.s_index().end();
        for (;it!=it_f;++it)
        {
          term.s_cf()=cf_type(it->g_cf());
          term.s_trig()=it->g_trig();
          term.s_flavour()=it->g_flavour();
          ancestor::insert(term);
        }
      }
/// Constructor from real and imaginary series.
      explicit complex(const real_ps &p, const real_ps &q):ancestor::base_pseries()
      {
        build_from_components(p,q);
      }
/// Constructor from real and imaginary series from filenames.
      explicit complex(const std::string &file1, const std::string &file2):ancestor::base_pseries()
      {
        build_from_components(real_ps(file1),real_ps(file2));
      }
  };
}
#endif
