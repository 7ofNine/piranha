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

#ifndef PIRANHA_POLYNOMIAL_H
#define PIRANHA_POLYNOMIAL_H

#include "base_polynomial.h"

namespace piranha
{
  template <class T>
    class polynomial : public base_polynomial<T,polynomial<T> >
  {
    public:
      typedef base_polynomial<T,polynomial> ancestor;
      typedef typename ancestor::m_type m_type;
      // Start INTERFACE definition for the real version.
      //-------------------------------------------------------
      /// Alias for the evaluation type.
      typedef typename m_type::eval_type eval_type;
      // Ctors.
      /// Default constructor.
      explicit polynomial():ancestor::base_polynomial()
        {}
      /// Constructor from int.
      explicit polynomial(int n):ancestor::base_polynomial(n)
        {}
      /// Constructor from double.
      explicit polynomial(const double &x):ancestor::base_polynomial(x)
        {}
      /// Constructor from psymbol
      explicit polynomial(const psymbol &psym):ancestor::base_polynomial(psym)
        {}
      /// Copy constructor.
      polynomial(const polynomial &p):ancestor::base_polynomial(p)
        {}
      explicit polynomial(const std::string &s):ancestor::base_polynomial(s)
        {}
      /// Destructor.
      ~polynomial()
        {}
      /// Assignment operator.
      polynomial &operator=(const polynomial &p)
      {
        ancestor::basic_assignment(p);
        return *this;
      }
      polynomial &operator*=(int n)
      {
        ancestor::mult_by_int(n);
        return *this;
      }
      polynomial &operator*=(const double &x)
      {
        ancestor::mult_by_double(x);
        return *this;
      }
      polynomial operator*(int n) const
      {
        polynomial retval(*this);
        retval*=n;
        return retval;
      }
      polynomial &operator+=(const polynomial &p)
      {
        ancestor::base_merge(p,true);
        return *this;
      }
      polynomial &operator-=(const polynomial &p)
      {
        ancestor::base_merge(p,false);
        return *this;
      }
      polynomial &operator/=(int n)
      {
        ancestor::ll_generic_integer_division(n);
        return *this;
      }
      polynomial operator-()
      {
        return (polynomial(*this)*=-1);
      }
      // Maths
      /// Bessel function of the first kind.
      polynomial besselJ(int n, const vector_psym_p &v) const
      {
        const unsigned int iterations=math::besselJ_series_limit(n,ancestor::t_eval(0,v));
        return math::pow_besselJ<polynomial,mpz_class>(n,*this,iterations);
      }
      polynomial pow(const double &x) const
      {
        polynomial retval(*this);
        retval.basic_pow(x);
        return retval;
      }
      // End INTERFACE definition.
      //-------------------------------------------------------
      polynomial &operator/=(const mpz_class &n)
      {
        ancestor::ll_generic_integer_division(n);
        return *this;
      }
  };
}


namespace std
{
  // COMPLEX COUNTERPART
  template <class T>
    struct complex<piranha::polynomial<T> > :
  public piranha::base_polynomial<complex<T>,complex<piranha::polynomial<T> > >
  {
    private:
      typedef complex self;
      typedef piranha::polynomial<T> real_type;
      typedef piranha::base_polynomial<complex<T>,complex<piranha::polynomial<T> > > ancestor;
      typedef typename ancestor::m_type m_type;
      typedef typename real_type::m_type real_m_type;
      typedef typename m_type::numerical_type numerical_type;
      typedef typename real_m_type::numerical_type real_numerical_type;
      // Access to iterators and indices of base.
      //typedef typename ancestor::expo_index expo_index;
      typedef typename ancestor::iterator iterator;
      typedef typename real_type::iterator real_iterator;
      // Start INTERFACE definition for the complex specialization.
      //-------------------------------------------------------
    public:
      typedef piranha::complex_double eval_type;
      // Ctors and dtor.
      explicit complex():ancestor::base_polynomial()
        {}
      explicit complex(const std::string &s):ancestor::base_polynomial(s)
        {}
      /// Constructor from real counterparts.
      /**
       * Takes monomials from r and i, turns them into complex monomials and inserts into self.
       */
      explicit complex(const real_type &r, const real_type &i)
      {
        // We do not need to care about widths because insertion in base_polynomial can cope with
        // larger and smaller polynomials.
        // In any case, this is a problem pseries should take care of, so maybe:
        // FIXME: assert equality between widths once we are using polynomials exclusively in pseries.
        real_iterator it, it_f;
        it_f=r.e_index().end();
        // Insert real part.
        for (it=r.e_index().begin();it!=it_f;++it)
        {
          ancestor::insert(m_type(*it));
        }
        // Insert imaginary part.
        it_f=i.e_index().end();
        m_type tmp;
        for (it=i.e_index().begin();it!=it_f;++it)
        {
          tmp=*it;
          tmp.numerical_cf()=numerical_type(real_numerical_type(0.),it->numerical_cf());
          ancestor::insert(tmp);
        }
      }
      bool operator==(const complex &p2) const
      {
        return ancestor::basic_comparison(p2);
      }
      explicit complex(const real_type &r)
      {
        // FIXME: share with ctor from real+imaginary?
        real_iterator it, it_f;
        it_f=r.e_index().end();
        // Insert real part.
        for (it=r.e_index().begin();it!=it_f;++it)
        {
          ancestor::insert(m_type(*it));
        }
      }
      explicit complex(const piranha::complex_double &c)
      {
        ancestor::insert(m_type(numerical_type(c)));
      }
      explicit complex(int n)
      {
        ancestor::insert(m_type(n));
      }
      // FIXME: where is this used?
      explicit complex(int n1, int n2)
      {
        m_type tmp(n1);
        ancestor::insert(tmp);
        tmp=m_type(n2);
        tmp.numerical_cf()=numerical_type(real_numerical_type(0.),tmp.numerical_cf().real());
        ancestor::insert(tmp);
      }
      /*explicit complex(const double &x):
          ancestor::simple_container(complex_double(x,0.))
      {}
      explicit complex(const double &x1, const double x2):
          ancestor::simple_container(complex_double(x1,x2))
      {}*/
      ~complex()
        {}
      real_type real() const
      {
        // FIXME: share this code with above ctors?
        // FIXME: use hinted insertion for improved performance.
        real_type retval;
        iterator it, it_f;
        it_f=ancestor::e_index().end();
        // Insert real part.
        for (it=ancestor::e_index().begin();it!=it_f;++it)
        {
          real_m_type tmp(it->numerical_cf().real(),it->rational_cf(),it->container());
          retval.insert(tmp);
        }
        return retval;
      }

      real_type imag() const
      {
        // FIXME: share this code with above ctors?
        // FIXME: share also with real().
        real_type retval;
        iterator it, it_f;
        it_f=ancestor::e_index().end();
        // Insert real part.
        for (it=ancestor::e_index().begin();it!=it_f;++it)
        {
          real_m_type tmp(it->numerical_cf().imag(),it->rational_cf(),it->container());
          retval.insert(tmp);
        }
        return retval;
      }
      // Setters.
      /// Like constructor from real.
      // FIXME: abstract and share.
      void set_real(const real_type &r)
      {
        ancestor::clear();
        real_iterator it, it_f;
        it_f=r.e_index().end();
        // Insert real part.
        for (it=r.e_index().begin();it!=it_f;++it)
        {
          ancestor::insert(m_type(*it));
        }
      }
      // FIXME: complete set/add API here.
      void set_imag(const real_type &i)
      {
        ancestor::clear();
        real_iterator it, it_f;
        it_f=i.e_index().end();
        // Insert imag part.
        for (it=i.e_index().begin();it!=it_f;++it)
        {
          ancestor::insert(m_type(numerical_type(real_numerical_type(0.),it->numerical_cf()),
            it->rational_cf(),it->container()));
        }
      }
      complex &operator=(const complex &c)
      {
        ancestor::basic_assignment(c);
        return *this;
      }
      complex operator-()
      {
        return (complex(*this)*=-1);
      }
      complex &operator/=(int n)
      {
        ancestor::ll_generic_integer_division(n);
        return *this;
      }
      complex &operator*=(int n)
      {
        ancestor::mult_by_int(n);
        return *this;
      }
      complex &operator*=(const double x)
      {
        ancestor::mult_by_double(x);
        return *this;
      }
      complex operator*(int n)
      {
        complex retval(*this);
        retval*=n;
        return retval;
      }
      complex &operator+=(const self &c)
      {
        ancestor::base_merge(c,true);
        return *this;
      }
      complex &operator-=(const complex &c)
      {
        ancestor::base_merge(c,false);
        return *this;
      }
      complex operator+(const complex &c) const
      {
        return (self(*this)+=c);
      }
      complex operator-(const complex &c) const
      {
        return (self(*this)-=c);
      }
  };
}
#endif
