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
#include "complex_toolbox.h"
#include "differential_toolbox.h"
#include "trigonometric_toolbox.h"

#define __COMMON_PS_TYPEDEFS(Cf,Trig,I)\
public:\
    typedef piranha::base_pseries<Cf,Trig,I> ancestor;\
    typedef typename ancestor::term_type term_type;\
    typedef typename ancestor::it_s_index it_s_index;\
    typedef typename ancestor::r_it_s_index r_it_s_index;\
    typedef typename ancestor::cf_type cf_type;

#define __COMMON_PS()                                                                             \
public:                                                                                           \
    ps():ancestor::base_pseries() {}                                                              \
    ps(const ps &p):ancestor::base_pseries(p) {}                                                  \
    explicit ps(const std::string &filename):ancestor::base_pseries(filename) {}                  \
    explicit ps(const cf_type &c):ancestor::base_pseries(c) {}                                    \
    ps copy() const {return ps(*this);}                                                           \


namespace piranha
{
  /// Default derived class for Poisson series.
  /**
   * This derived class provides support for complex arithmetic and for trigonometric operations on
   * series.
   */
  template <class Cf, class Trig, template <class,class> class I>
  class ps:
        public base_pseries<Cf,Trig,I>,
        public common_trig_toolbox<ps<Cf,Trig,I>,ps<Cf,Trig,I> >,
        public real_trig_toolbox<ps<Cf,Trig,I> >,
        public differential_toolbox<ps<Cf,Trig,I> >
    {
      //|-----------------------|
      //|Typedef Specializations|
      //|-----------------------|
      __COMMON_PS_TYPEDEFS(Cf,Trig,I);
    public:
      /// Alias for self.
      typedef ps real_ps;
      /// Alias for complex coefficient type.
      typedef std::complex<cf_type> complex_cf_type;
      /// Alias for complex counterpart.
      typedef std::complex<ps<cf_type,Trig,I> > complex_ps;
      // Common stuff.
      __COMMON_PS();
      __PS_OPERATORS(ps,ancestor);
      // FIXME: cram this into operator toolbox.
      ps &operator/=(int n)
      {
        ancestor::basic_div_by_int(n);
        return *this;
      }
      //|----------------------|
      //|Public Specializations|
      //|----------------------|
    public:
      /// Constructor from int.
      explicit ps(int n)
      {
        term_type tmp=term_type(cf_type(n));
        ancestor::insert(tmp);
      }
      /// Constructor from double.
      explicit ps(const double &x)
      {
        term_type tmp=term_type(cf_type(x));
        ancestor::insert(tmp);
      }
      /// Constructor from psymbol.
      explicit ps(const psymbol &psym):ancestor::base_pseries(psym)
      {}
    };
}


namespace std
  {
  // COMPLEX COUNTERPART
  /// Complex specialization for default derived class.
  template <class Cf,class Trig,template <class,class> class I>
  struct complex<piranha::ps<Cf,Trig,I> >:
        public piranha::base_pseries <complex<Cf>,Trig,I>,
        public piranha::common_trig_toolbox<complex<piranha::ps<Cf,Trig,I> >,piranha::ps<Cf,Trig,I> >,
        public piranha::complex_toolbox<piranha::ps<Cf,Trig,I> >,
        public differential_toolbox<complex<piranha::ps<Cf,Trig,I> > >
    {
public:
      /// Alias for the ancestor.
      typedef piranha::base_pseries<complex<Cf>,Trig,I> ancestor;
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
      typedef piranha::ps<Cf,Trig,I> real_ps;
      /// Alias for real coefficient.
      typedef Cf real_cf_type;
      // Common stuff
      // ---------------------------------------------------------
      // FIXME: this is shared with above, find a way to share.
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
      /// Copy to another complex series.
      complex copy() const
        {
          return complex(*this);
        }
      // ---------------------------------------------------------
      //__PS_OPERATORS(ps,ancestor);
      //|----------------------|
      //|Public Specializations|
      //|----------------------|
public:
      /// Constructor from pair of real coefficients.
      explicit complex(const real_cf_type &a, const real_cf_type &b)
      {
        term_type tmp(cf_type(a,b));
        ancestor::insert(tmp);
      }
      /// Constructor from complex.
      explicit complex(const complex_double &c)
      {
        term_type tmp=term_type(cf_type(c));
        ancestor::insert(tmp);
      }
      /// Constructor from pair of doubles.
      explicit complex(const double &a, const double &b)
      {
        term_type tmp(cf_type(a,b));
        ancestor::insert(tmp);
      }
      /// Constructor from real series.
      // WARNING: here and below we are discarding lin_args.
      explicit complex(const real_ps &p):ancestor::base_pseries()
      {
        p_assert(ancestor::merge_args(p));
        term_type term;
        typename real_ps::it_s_index it=p.s_index().begin(), it_f=p.s_index().end();
        for (;it!=it_f;++it)
          {
            term.c()=cf_type(it->c());
            term.trig_args()=it->trig_args();
            term.flavour()=it->flavour();
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
      // OPERATORS
      // FIXME: cram this into operator toolbox.
      complex &operator/=(int n)
      {
        ancestor::basic_div_by_int(n);
        return *this;
      }
      template <class Cf2,class Trig2,template <class,class> class I2>
      complex &operator*=(const complex<piranha::ps<Cf2,Trig2,I2> > &p)
      {
        ancestor::basic_ps_mult(p);
        return *this;
      }
      template <class Cf2,class Trig2,template <class,class> class I2>
      complex &operator*=(const piranha::ps<Cf2,Trig2,I2> &p)
      {
        ancestor::basic_ps_mult(p);
        return *this;
      }
      complex &operator*=(int n)
      {
        ancestor::mult_by_int(n);
        return *this;
      }
      complex &operator*=(const double &x)
      {
        ancestor::generic_mult(x);
        return *this;
      }
      complex &operator*=(const long double &x)
      {
        ancestor::generic_mult(x);
        return *this;
      }
      template <class Cf2,class Trig2,template <class,class> class I2>
      complex &operator+=(const complex<piranha::ps<Cf2,Trig2,I2> > &p)
      {
        ancestor::merge_with(p);
        return *this;
      }
      template <class Cf2,class Trig2,template <class,class> class I2>
      complex &operator+=(const piranha::ps<Cf2,Trig2,I2> &p)
      {
        ancestor::merge_with(p);
        return *this;
      }
      complex &operator+=(int x)
      {
        ancestor::generic_merge(x);
        return *this;
      }
      complex &operator+=(const double &x)
      {
        ancestor::generic_merge(x);
        return *this;
      }
      complex &operator+=(const long double &x)
      {
        ancestor::generic_merge(x);
        return *this;
      }
      template <class Cf2,class Trig2,template <class,class> class I2>
      complex operator+(const complex<piranha::ps<Cf2,Trig2,I2> > &p) const
        {
          complex retval(*this);
          retval+=p;
          return retval;
        }
      template <class Cf2,class Trig2,template <class,class> class I2>
      complex operator+(const piranha::ps<Cf2,Trig2,I2> &p) const
        {
          complex retval(*this);
          retval+=p;
          return retval;
        }
      complex operator+(const double &x) const
        {
          complex retval(*this);
          retval+=x;
          return retval;
        }
    };
}
#endif
