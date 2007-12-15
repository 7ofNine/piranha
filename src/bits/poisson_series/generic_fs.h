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

#include "base_classes/base_pseries.h"
#include "toolboxes/operators_toolbox.h"
#include "toolboxes/fourier_multiplication_toolbox.h"
#include "toolboxes/math_toolbox.h"
#include "toolboxes/complex_toolbox.h"
#include "toolboxes/differential_toolbox.h"
#include "toolboxes/trigonometric_toolbox.h"

namespace piranha
{
/// Derived class for Fourier series.
/**
 * This derived class provides support for complex arithmetic and for trigonometric operations on
 * Fourier series. Multiplication is norm-truncated.
 */
  template <class Cf, class Trig, template <class,class> class Term,
    template <class,class, template <class, class> class > class I, class Allocator=std::allocator<char> >
    class generic_fs:
    public base_pseries<Cf,Trig,Term,I,generic_fs<Cf,Trig,Term,I,Allocator>,Allocator>,
    public real_operators_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public math_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public trigonometric_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public differential_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >,
    public fourier_multiplication_toolbox<generic_fs<Cf,Trig,Term,I,Allocator> >
  {
      template <class Derived> friend class real_operators_toolbox;
      template <class Derived> friend class complex_toolbox;
      template <class Derived> friend class base_trigonometric_toolbox;
      template <class Derived> friend class fourier_multiplication_toolbox;
      template <class Derived> friend class math_toolbox;
      template <class Derived> friend class differential_toolbox;
    public:
/// Alias for parent class.
      typedef piranha::base_pseries<Cf, Trig, Term, I, generic_fs<Cf,Trig,Term,I,Allocator>,Allocator> ancestor;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Empty constructor.
      generic_fs():ancestor::base_pseries() {}
/// Copy constructor.
      generic_fs(const generic_fs &p):ancestor::base_pseries(p) {}
/// Constructor from file.
      explicit generic_fs(const std::string &filename):ancestor::base_pseries(filename) {}
/// Constructor from int.
      explicit generic_fs(int n):ancestor::base_pseries(n) {}
/// Constructor from double.
      explicit generic_fs(const double &x):ancestor::base_pseries(x) {}
/// Constructor from psymbol.
      explicit generic_fs(const psymbol &psym, psymbol::type ptype):ancestor::base_pseries(psym,ptype) {}
/// Constructor from coefficient and model.
      explicit generic_fs(const cf_type &c, const generic_fs &model):ancestor::base_pseries(c,model) {}
  };
}

namespace std
{
/// Complex specialization for Fourier series derived class.
  template <class Cf, class Trig, template <class,class> class Term,
    template <class,class, template <class, class> class > class I, class Allocator>
    struct complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >:
    public piranha::base_pseries<complex<Cf>,Trig,Term,I,complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,Allocator>,
    public piranha::complex_operators_toolbox<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,
    public piranha::trigonometric_toolbox<complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > >,
    public piranha::complex_toolbox<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,
    public piranha::differential_toolbox<complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > >,
    public piranha::fourier_multiplication_toolbox<complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > >
  {
      template <class Derived> friend class piranha::fourier_multiplication_toolbox;
      template <class Derived> friend class piranha::base_trigonometric_toolbox;
      template <class Derived> friend class piranha::trigonometric_toolbox;
      template <class Derived> friend class piranha::complex_toolbox;
      template <class Derived> friend class differential_toolbox;
    public:
/// Alias for ancestor.
      typedef piranha::base_pseries<complex<Cf>,Trig,Term,I,complex<piranha::generic_fs<Cf,Trig,Term,I,Allocator> >,
        Allocator> ancestor;
/// Alias for coefficient type.
      typedef typename ancestor::cf_type cf_type;
/// Alias for real counterpart.
      typedef piranha::generic_fs<Cf,Trig,Term,I,Allocator> value_type;
/// Alias for real coefficient.
      typedef Cf real_cf_type;
/// Alias for complex toolbox.
      typedef piranha::complex_toolbox<piranha::generic_fs<Cf,Trig,Term,I,Allocator> > complex_toolbox;
// Base ctors.
/// Default constructor.
      explicit complex():ancestor::base_pseries() {}
/// Copy constructor.
      complex(const complex &p):ancestor::base_pseries(p),complex_toolbox::complex_toolbox(p) {}
/// Constructor from filename.
      explicit complex(const std::string &filename):ancestor::base_pseries(filename) {}
/// Constructor from int.
      explicit complex(int n):ancestor::base_pseries(n) {}
/// Constructor from double.
      explicit complex(const double &x):ancestor::base_pseries(x) {}
/// Constructor from psymbol.
      explicit complex(const piranha::psymbol &psym, piranha::psymbol::type ptype):ancestor::base_pseries(psym,ptype) {}
/// Constructor from coefficient and model.
      explicit complex(const cf_type &c, const complex &model):ancestor::base_pseries(c,model) {}
// Complex specific ctors.
/// Constructor from complex integer.
      explicit complex(const complex<int> &n):ancestor::base_pseries(),complex_toolbox(n) {}
/// Constructor from real and imaginary integers.
      explicit complex(int r, int i):ancestor::base_pseries(),complex_toolbox(r,i) {}
/// Constructor from complex double.
      explicit complex(const complex<double> &x):ancestor::base_pseries(),complex_toolbox(x) {}
/// Constructor from real and imaginary doubles.
      explicit complex(const double &r, const double &i):ancestor::base_pseries(),complex_toolbox(r,i) {}
/// Constructor from value_type.
      explicit complex(const value_type &p):ancestor::base_pseries(),complex_toolbox(p) {}
/// Constructor from real and imaginary series.
      explicit complex(const value_type &p, const value_type &q):ancestor::base_pseries(),complex_toolbox(p,q) {}
  };
}
#endif
