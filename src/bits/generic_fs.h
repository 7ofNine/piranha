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

#ifndef PIRANHA_GENERIC_FS_H
#define PIRANHA_GENERIC_FS_H

#include <complex>

#include "../base_classes/base_pseries.h"
#include "../toolboxes/operators_toolbox.h"
#include "../toolboxes/fourier_multiplication_toolbox.h"
#include "../toolboxes/math_toolbox.h"
#include "../toolboxes/complex_toolbox.h"
#include "../toolboxes/differential_toolbox.h"
#include "../toolboxes/trigonometric_toolbox.h"

namespace piranha
{
/// Derived class for Fourier series.
/**
 * This derived class provides support for complex arithmetic and for trigonometric operations on
 * Fourier series. Multiplication is norm-truncated.
 */
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I, class Allocator=std::allocator<char> >
    class generic_fs:
    public base_pseries<Cf,Trig,Term,I,generic_fs<Cf,Trig,Term,I,Allocator>,Allocator>,
    public real_operators_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public math_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public trigonometric_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public differential_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public fourier_multiplication_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >
  {
    public:
/// Alias for parent class.
      typedef piranha::base_pseries<Cf, Trig, Term, I, generic_fs<Cf,Trig,Term,I,Allocator>,Allocator> ancestor;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Alias for term type.
      typedef typename ancestor::term_type term_type;
/// Alias for self.
      typedef generic_fs<Cf,Trig,Term,I,Allocator> self;
// Ctors.
      generic_fs():ancestor::base_pseries()
        {}
      generic_fs(const generic_fs &p):ancestor::base_pseries(p)
        {}
      explicit generic_fs(const std::string &filename):ancestor::base_pseries(filename)
        {}
      explicit generic_fs(const cf_type &c, const generic_fs &model):ancestor::base_pseries(c,model)
        {}
/// Constructor from int.
      explicit generic_fs(int n)
      {
        ancestor::generic_builder(n);
      }
/// Constructor from double.
      explicit generic_fs(const double &x)
      {
        ancestor::generic_builder(x);
      }
/// Constructor from psymbol.
      explicit generic_fs(const psymbol &psym, psymbol::type ptype):ancestor::base_pseries(psym,ptype)
        {}
  };
}


namespace std
{
// COMPLEX COUNTERPART
/// Complex specialization for Fourier series derived class.
  template <class Cf, class Trig, template <class,class> class Term, template <class,class, template <class, class> class > class I, class Allocator>
    struct complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >:
    public piranha::base_pseries<complex<Cf>,Trig,Term,I,complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,Allocator>,
    public piranha::complex_operators_toolbox<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,
    public piranha::trigonometric_toolbox<complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > >,
    public piranha::complex_toolbox<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,
    public piranha::differential_toolbox<complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > >,
    public piranha::fourier_multiplication_toolbox<complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > >
  {
    public:
/// Alias for the ancestor.
      typedef piranha::base_pseries<complex<Cf>,Trig,Term,I,complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,
        Allocator> ancestor;
/// Alias for term type.
      typedef typename ancestor::term_type term_type;
/// Alias for the iterator to sorted index.
      typedef typename ancestor::it_s_index it_s_index;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Alias for real counterpart.
      typedef piranha::generic_fs<Cf,Trig,Term,I,Allocator> real_type;
/// Alias for real coefficient.
      typedef Cf real_cf_type;
/// Alias for self.
      typedef complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > self;
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
      explicit complex(const cf_type &c, const complex &model):ancestor::base_pseries(c,model)
        {}
/// Constructor from pair of real coefficients.
/*explicit complex(const real_cf_type &a, const real_cf_type &b):
  ancestor::base_pseries(cf_type(a,b))
  {}*/
/// Constructor from complex.
      explicit complex(const piranha::complex_double &c)
      {
        ancestor::generic_builder(c);
      }
/// Constructor from pair of doubles.
      explicit complex(const double &a, const double &b)
      {
        ancestor::generic_builder(piranha::complex_double(a,b));
      }
/// Constructor from real series.
// FIXME: here and below we are discarding lin_args.
// TODO: can we re-use some function from complex_toolbox to achieve this result?
// If so, ditch the term_type typedef which is used just here. Also the iterator typedef..
      explicit complex(const real_type &p)
      {
        action_assert(ancestor::merge_args(p));
        term_type term;
        typename real_type::it_s_index it=p.g_s_index().begin(), it_f=p.g_s_index().end();
        for (;it!=it_f;++it)
        {
          *term.s_cf()=cf_type(*it->g_cf());
          *term.s_trig()=*it->g_trig();
          term.s_flavour()=it->g_flavour();
          ancestor::insert(term);
        }
      }
/// Constructor from real and imaginary series.
      explicit complex(const real_type &p, const real_type &q)
      {
        build_from_components(p,q);
      }
/// Constructor from real and imaginary series from filenames.
      explicit complex(const std::string &file1, const std::string &file2)
      {
        build_from_components(real_type(file1),real_type(file2));
      }
  };
}
#endif
