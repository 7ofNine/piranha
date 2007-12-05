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
#include "../bits/concepts/basic_pseries_coefficient_concept.h"

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
    public basic_pseries_coefficient_concept<double_cf>,
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
      using ancestor::add;
      using ancestor::subtract;
      using ancestor::mult_by;
      using ancestor::mult_by_self;
      using ancestor::divide_by;
// Start implementation of basic pseries coefficient interface.
//------------
// Ctors and dtor.
/// Empty constructor.
      explicit double_cf():ancestor::numerical_container() {}
/// Constructor from string.
      explicit double_cf(const std::string &s):ancestor::numerical_container(s) {}
/// Constructor from symbol.
      explicit double_cf(const psymbol &):ancestor::numerical_container()
      {
        std::cout << "WARNING: building numerical coefficient from psymbol." << std::endl;
      }
/// Constructor from integer.
      explicit double_cf(int val):ancestor::numerical_container(val) {}
/// Constructor from double.
      explicit double_cf(const double &val):ancestor::numerical_container(val) {}
/// Copy constructor.
      double_cf(const self &dc):ancestor::numerical_container(dc) {}
/// Destructor.
      ~double_cf() {}
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
      self pow(const double &y) const
      {
        self retval;
        retval.s_value()=std::pow(g_value(),y);
        return retval;
      }
// Needed operators.
      self &operator=(const self &val2)
      {
        s_value()=val2.g_value();
        return *this;
      }
// End implementation of basic pseries coefficient interface.
//------------
// Start implementation of trigonometric pseries coefficient interface.
// Used in:
// - trigonometric toolbox,
//------------
// Maths
/// Bessel function of the first kind.
/**
 * Uses C standard library call.
 */
      self besselJ(int n, const vector_psym_p &) const
      {
        self retval;
        retval.s_value()=math::besselJ(n,g_value());
        return retval;
      }
// End implementation of trigonometric pseries coefficient interface.
//------------
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
    public piranha::complex_basic_pseries_coefficient_concept<piranha::double_cf>,
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
      using ancestor::add;
      using ancestor::subtract;
      using ancestor::mult_by;
      using ancestor::mult_by_self;
      using ancestor::divide_by;
// Start implementation of basic pseries coefficient interface.
//------------
// Ctors and dtor.
      explicit complex():ancestor::numerical_container() {}
      explicit complex(const std::string &s):ancestor::numerical_container(s) {}
      explicit complex(const piranha::psymbol &):ancestor::numerical_container() {}
      explicit complex(int n):ancestor::numerical_container(n) {}
      explicit complex(const double &x):ancestor::numerical_container(x) {}
// TODO: put those ctors in toolbox and use from there.
/// Constructor from pair of ints.
      explicit complex(int r, int i):ancestor::numerical_container()
      {
        s_value().real()=r;
        s_value().imag()=i;
      }
/// Constructor from complex integer.
      explicit complex(const std::complex<int> &c):ancestor::numerical_container()
      {
        s_value()=c;
      }
/// Constructor from pair of doubles.
      explicit complex(const double &r, const double &i):ancestor::numerical_container()
      {
        s_value().real()=r;
        s_value().imag()=i;
      }
/// Constructor from complex double.
      explicit complex(const std::complex<double> &c):ancestor::numerical_container()
      {
        s_value()=c;
      }
/// Constructor from real type.
      explicit complex(const value_type &r):ancestor::numerical_container()
      {
        s_value().real()=r.g_value();
      }
/// Constructor from real and imaginary parts.
      explicit complex(const value_type &r, const value_type &i):ancestor::numerical_container()
      {
        s_value().real()=r.g_value();
        s_value().imag()=i.g_value();
      }
/// Copy constructor.
      complex(const complex &c):ancestor::numerical_container(c) {}
      ~complex() {}
// Getters.
      double norm(const piranha::vector_psym_p &) const
      {
        return abs();
      }
      value_type real() const
      {
        value_type retval;
        retval.s_value()=g_value().real();
        return retval;
      }
      value_type imag() const
      {
        value_type retval;
        retval.s_value()=g_value().imag();
        return retval;
      }
// Setters.
/// Set value_ to be a real only value.
      void set_real(const value_type &r)
      {
        s_value()=r.g_value();
      }
/// Set value_ to be an imag only value.
      void set_imag(const value_type &i)
      {
        s_value()=piranha::complex_double(0,i.g_value());
      }
// Probing.
      bool is_zero(const piranha::vector_psym_p &) const
      {
        return (abs()<piranha::settings_manager::numerical_zero());
      }
// Maths.
      template <class DerivedPs>
        self &mult_by_self(const value_type &x, const DerivedPs &)
      {
        return mult_by_generic(x.g_value());
      }
      self &mult_by(const std::complex<int> &c)
      {
        return mult_by_generic(c);
      }
      self &mult_by(const std::complex<double> &c)
      {
        return mult_by_generic(c);
      }
// Operators.
      self &operator=(const self &val2)
      {
        s_value()=val2.g_value();
        return *this;
      }
      self &operator=(const value_type &r2)
      {
        s_value()=r2.g_value();
        return *this;
      }
// End implementation of complex basic pseries coefficient interface.
//------------
  };

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
