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

#ifndef PIRANHA_DOUBLE_CF_H
#define PIRANHA_DOUBLE_CF_H

#include <complex>

#include "../bits/common_typedefs.h"
#include "../bits/math.h"                                 // besselJ.
#include "../base_classes/numerical_container.h"
#include "../bits/concepts/pseries_coefficient_concept.h"

namespace piranha
{
/// Double-precision numerical coefficient.
/**
 * This class can be used as coefficient in Poisson series. It encapsulate a double precision
 * numerical value and provides the means to access it.
 *
 * A set of operators is provided to enable interoperability with basic numerical data types.
 */
  class double_cf :
    public pseries_coefficient_concept<double_cf>,
    public numerical_container<double,double_cf>
  {
    public:
/// Alias for self.
      typedef double_cf self;
/// Alias for the parent class.
      typedef numerical_container<double,double_cf> ancestor;
// This is necessary because some moethods are present in concept _and_ in numerical container.
// We avoid the dreaded diamond problem by explicitly stating which functions to use.
      using ancestor::swap;
      using ancestor::print_plain;
      using ancestor::print_latex;
      using ancestor::is_compatible;
      using ancestor::checkup;
      using ancestor::invert_sign;
      using ancestor::t_eval;
      using ancestor::add_self;
      using ancestor::subtract_self;
      using ancestor::mult_by_int;
      using ancestor::mult_by_double;
      using ancestor::mult_by_generic;
      using ancestor::mult_by_self;
      using ancestor::divide_by_int;
      using ancestor::divide_by_generic;
// Start INTERFACE definition for the real version.
//-------------------------------------------------------
// Ctors and dtor.
/// Empty constructor.
      explicit double_cf():ancestor::numerical_container()
        {}
/// Constructor from double.
      explicit double_cf(const double &val):ancestor::numerical_container(val)
        {}
/// Constructor from integer.
      explicit double_cf(int val):ancestor::numerical_container(val)
        {}
/// Constructor from string.
      explicit double_cf(const std::string &s):ancestor::numerical_container(s)
        {}
/// Constructor from symbol.
      explicit double_cf(const psymbol &):ancestor::numerical_container()
      {
        std::cout << "WARNING: building numerical coefficient from psymbol." << std::endl;
      }
/// Copy constructor.
      double_cf(const self &dc):ancestor::numerical_container(dc)
        {}
/// Destructor.
      ~double_cf()
        {}
// Getters
/// Calculate norm (absolute value).
      double norm(const vector_psym_p &) const
      {
        return abs();
      }
// Probing
/// Is value zero?
/**
 * If value is less than settings_manager::numerical_zero() in absolute value it is considered
 * to be zero.
 */
      bool is_zero(const vector_psym_p &) const
      {
        return (abs()<settings_manager::numerical_zero());
      }
// Maths
/// Bessel function of the first kind.
/**
 * Uses C standard library call.
 */
      self besselJ(int n, const vector_psym_p &) const
      {
        return self(math::besselJ(n,value_));
      }
      self pow(const double &y) const
      {
        return self(std::pow(value_,y));
      }
// Needed operators.
      self &operator=(const self &val2)
      {
        ancestor::value_=val2.value_;
        return *this;
      }
// End INTERFACE definition.
//-------------------------------------------------------
// Interface for monomial.
      double abs() const
      {
        return std::abs(value_);
      }
/// Test whether two double_cf are equal or opposite in sign.
      bool equal_or_opposite(const self &val2) const
      {
        return (std::abs(abs()-val2.abs())<settings_manager::numerical_zero());
      }
/// Test whether two double_cf have the same sign.
      bool same_sign(const self &val2) const
      {
        return ((value_>0)==(val2.value()>0));
      }
/// Test whether value is unity.
      bool is_unity() const
      {
        return (std::abs(abs()-1)<settings_manager::numerical_zero());
      }
/// Test whether value is negative.
      bool is_negative() const
      {
        return (value_<0);
      }
/// Test whether value is positive.
      bool is_positive() const
      {
        return (value_>0);
      }
/// Calculate norm (absolute value).
      double norm() const
      {
        return abs();
      }
      self inv() const
      {
        return self(1./value_);
      }
      bool operator<(const double &x) const
      {
        return (value_<x);
      }
// ancestor::value() is part of the interface too.
//-------------------------------------------------------
  };

  inline std::istream &operator>>(std::istream &is, double_cf &dc)
  {
    std::string tmp;
    getline(is,tmp);
    dc.value()=utils::lexical_converter<double>(tmp);
    return is;
  }

  inline std::ostream &operator<<(std::ostream &os, const double_cf &dc)
  {
    os << dc.value();
    return os;
  }
}


namespace std
{
  template <>
    struct complex<piranha::double_cf>:
    public complex<piranha::pseries_coefficient_concept<piranha::double_cf> >,
    public piranha::numerical_container<piranha::complex_double,complex<piranha::double_cf> >
  {
      typedef piranha::numerical_container<piranha::complex_double,complex<piranha::double_cf> > ancestor;
      using ancestor::swap;
      using ancestor::print_plain;
      using ancestor::print_latex;
      using ancestor::is_compatible;
      using ancestor::checkup;
      using ancestor::invert_sign;
      using ancestor::t_eval;
      using ancestor::add_self;
      using ancestor::subtract_self;
      using ancestor::mult_by_int;
      using ancestor::mult_by_double;
      using ancestor::mult_by_generic;
      using ancestor::mult_by_self;
      using ancestor::divide_by_int;
      using ancestor::divide_by_generic;
    public:
      typedef complex self;
      typedef piranha::double_cf double_type;
      typedef double_type real_self;
// Start INTERFACE definition for the complex specialization. FIXME: is this different from
// the above???
//-------------------------------------------------------
// Ctors and dtor.
      explicit complex():ancestor::numerical_container()
        {}
      explicit complex(const std::string &s):ancestor::numerical_container(s)
        {}
      explicit complex(const double_type &r, const double_type &i):
      ancestor::numerical_container(piranha::complex_double(r.value(),i.value()))
        {}
      explicit complex(const double_type &r):
      ancestor::numerical_container(piranha::complex_double(r.value(),0.))
        {}
      explicit complex(const piranha::complex_double &c):ancestor::numerical_container(c)
        {}
      explicit complex(int n):ancestor::numerical_container(piranha::complex_double(n,0.))
        {}
      explicit complex(int n1, int n2):ancestor::numerical_container(piranha::complex_double(n1,n2))
        {}
      explicit complex(const double &x):
      ancestor::numerical_container(piranha::complex_double(x,0.))
        {}
      explicit complex(const double &x1, const double x2):
      ancestor::numerical_container(piranha::complex_double(x1,x2))
        {}
      ~complex()
        {}
// Getters.
      double norm(const piranha::vector_psym_p &) const
      {
        return abs();
      }
      double_type real() const
      {
        return double_type(value_.real());
      }
      double_type imag() const
      {
        return double_type(value_.imag());
      }
// Setters.
/// Set value_ to be a real only value.
      void set_real(const double_type &r)
      {
        value_=r.value();
      }
/// Set value_ to be an imag only value.
      void set_imag(const double_type &i)
      {
        value_=piranha::complex_double(0,i.value());
      }
// Probing.
      bool is_zero(const piranha::vector_psym_p &) const
      {
        return (abs()<piranha::settings_manager::numerical_zero());
      }
// Operators.
      self &operator=(const self &val2)
      {
        ancestor::value_=val2.value_;
        return *this;
      }
      self &operator=(const double_type &r2)
      {
        ancestor::value_=r2.value();
        return *this;
      }
// End INTERFACE definition.
//-------------------------------------------------------
// Interface for monomial.
/// Absolute value.
      double abs() const
      {
        return std::abs(value_);
      }
/// Test whether two complex double_cf are equal or opposite in sign.
      bool equal_or_opposite(const self &val2) const
      {
        return (std::abs(value_+val2.value())<piranha::settings_manager::numerical_zero() ||
          std::abs(value_-val2.value())<piranha::settings_manager::numerical_zero());
      }
/// Test whether two complex double_cf are in the same quadrant of the complex plane.
      bool same_sign(const self &val2) const
      {
        return ((value_.real()>=0 && val2.value().real()>=0) || (value_.real()<0 && val2.value().real()<0)) &&
          ((value_.imag()>=0 && val2.value().imag()>=0) || (value_.imag()<0 && val2.value().imag()<0));
      }
/// Test whether value is real unity.
      bool is_unity() const
      {
        return (is_real() &&
          std::abs(std::abs(value_.real())-1)<piranha::settings_manager::numerical_zero());
      }
/// Test whether value is negative.
      bool is_negative() const
      {
        return (value_.real()<0 && is_real());
      }
/// Test whether value is positive.
      bool is_positive() const
      {
        return (value_.real()>0 && is_real());
      }
// ancestor::value() is part of the interface too.
//-------------------------------------------------------
    private:
      bool is_real() const
      {
        return (std::abs(value_.imag())<piranha::settings_manager::numerical_zero());
      }
  }
  ;

// Overloads for I/O operators.
  inline istream &operator>>(istream &is, complex<piranha::double_cf> &dc)
  {
    string tmp;
    getline(is,tmp);
    dc.value()=piranha::utils::lexical_converter<piranha::complex_double>(tmp);
    return is;
  }

  inline ostream &operator<<(ostream &os, const complex<piranha::double_cf> &dc)
  {
    os << dc.value();
    return os;
  }
}
#endif
