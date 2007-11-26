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
/// Alias for self.
      typedef double_cf self;
/// Alias for the parent class.
      typedef numerical_container<double,double_cf> ancestor;
    public:
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
/// Constructor from string.
      explicit double_cf(const std::string &s):ancestor::numerical_container(s)
        {}
/// Constructor from double.
       explicit double_cf(const double &val):ancestor::numerical_container(val)
         {}
/// Constructor from integer.
       explicit double_cf(int val):ancestor::numerical_container(val)
         {}
/// Constructor from symbol.
      explicit double_cf(const psymbol &):ancestor::numerical_container()
      {
        std::cout << "WARNING: building numerical coefficient from psymbol." << std::endl;
      }
/// Generic constructor.
      template <class T>
        explicit double_cf(const T &x):ancestor::numerical_container(x)
      {}
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
        self retval;
        retval.value_=math::besselJ(n,value_);
        return retval;
      }
      self pow(const double &y) const
      {
        self retval;
        retval.value_=std::pow(value_,y);
        return retval;
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
        return ((value_>0)==(val2.g_value()>0));
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
        if (std::abs(value_)<settings_manager::numerical_zero())
        {
          std::cout << "Error: division by zero." << std::endl;
          std::abort();
        }
        self retval;
        retval.value_=1./value_;
        return retval;
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
    dc.s_value()=utils::lexical_converter<double>(tmp);
    return is;
  }

  inline std::ostream &operator<<(std::ostream &os, const double_cf &dc)
  {
    os << dc.g_value();
    return os;
  }
}

namespace std
{
  template <>
    struct complex<piranha::double_cf>:
    public piranha::complex_pseries_coefficient_concept<complex<piranha::double_cf> >,
    public piranha::numerical_container<piranha::complex_double,complex<piranha::double_cf> >
  {
    private:
      typedef piranha::numerical_container<piranha::complex_double,complex<piranha::double_cf> > ancestor;
      typedef complex self;
      typedef piranha::double_cf value_type;
    public:
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
// Start INTERFACE definition for the complex specialization. FIXME: is this different from
// the above???
//-------------------------------------------------------
// Ctors and dtor.
      explicit complex():ancestor::numerical_container()
      {}
      explicit complex(const std::string &s):ancestor::numerical_container(s)
      {}
//       explicit complex(const value_type &r, const value_type &i):
//       ancestor::numerical_container(piranha::complex_double(r.g_value(),i.g_value()))
//         {}
//       explicit complex(const value_type &r):ancestor::numerical_container(r)
//       {}
      explicit complex(const piranha::complex_double &c):ancestor::numerical_container(c)
      {}
      explicit complex(int n):ancestor::numerical_container(piranha::complex_double((double)n,0.))
      {}
//       explicit complex(int n1, int n2):ancestor::numerical_container(piranha::complex_double(n1,n2))
//       {}
      explicit complex(const double &x):ancestor::numerical_container(x)
      {}
/// Constructor from real type.
      explicit complex(const value_type &r)
      {
        ancestor::s_value().real()=r.g_value();
        ancestor::s_value().imag()=0;
      }
/// Constructor from real and imaginary parts.
      explicit complex(const value_type &r, const value_type &i)
      {
        ancestor::s_value().real()=r.g_value();
        ancestor::s_value().imag()=i.g_value();
      }
/// Generic constructor.
      template <class T>
        explicit complex(const T &x):ancestor::numerical_container(x)
      {}
/*      explicit complex(const double &x1, const double x2):
      ancestor::numerical_container(piranha::complex_double(x1,x2))
      {}*/
// FIXME: add copy ctor, assignment, etc.
/// Copy constructor.
      complex(const complex &c):ancestor::numerical_container(c)
      {}
      ~complex()
      {}
// Getters.
      double norm(const piranha::vector_psym_p &) const
      {
        return abs();
      }
      value_type real() const
      {
        value_type retval;
        retval.s_value()=value_.real();
        return retval;
      }
      value_type imag() const
      {
        value_type retval;
        retval.s_value()=value_.imag();
        return retval;
      }
// Setters.
/// Set value_ to be a real only value.
      void set_real(const value_type &r)
      {
        value_=r.g_value();
      }
/// Set value_ to be an imag only value.
      void set_imag(const value_type &i)
      {
        value_=piranha::complex_double(0,i.g_value());
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
      self &operator=(const value_type &r2)
      {
        ancestor::value_=r2.g_value();
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
        return (std::abs(value_+val2.g_value())<piranha::settings_manager::numerical_zero() ||
          std::abs(value_-val2.g_value())<piranha::settings_manager::numerical_zero());
      }
/// Test whether two complex double_cf are in the same quadrant of the complex plane.
      bool same_sign(const self &val2) const
      {
        return ((value_.real()>=0 && val2.g_value().real()>=0) || (value_.real()<0 && val2.g_value().real()<0)) &&
          ((value_.imag()>=0 && val2.g_value().imag()>=0) || (value_.imag()<0 && val2.g_value().imag()<0));
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
    dc.s_value()=piranha::utils::lexical_converter<piranha::complex_double>(tmp);
    return is;
  }

  inline ostream &operator<<(ostream &os, const complex<piranha::double_cf> &dc)
  {
    os << dc.g_value();
    return os;
  }
}
#endif
