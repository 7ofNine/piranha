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
  class double_cf:
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
      using ancestor::norm;
      using ancestor::is_zero;
      using ancestor::increase_size;
      using ancestor::append_args;
      using ancestor::prepend_args;
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
// Needed operators.
      self &operator=(const self &val2)
      {
        return assign_self(val2);
      }
// End implementation of basic pseries coefficient interface.
//------------
// Start implementation of trigonometric pseries coefficient interface.
// Used in:
// - trigonometric toolbox,
//------------
// TODO: Move into own toolbox and concept.
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
// Start implementation of power-enabled pseries coefficient interface.
// TODO: Move into own toolbox and concept.
      self pow(const double &y) const
      {
        self retval;
        retval.s_value()=std::pow(g_value(),y);
        return retval;
      }
  };
}

namespace std
{
  template <>
    struct complex<piranha::double_cf>:
    public piranha::complex_basic_pseries_coefficient_concept<piranha::double_cf>,
    public piranha::numerical_container<piranha::complex_double,complex<piranha::double_cf> >,
    public piranha::numerical_container_complex_toolbox<piranha::double_cf>
  {
    private:
      typedef piranha::numerical_container<piranha::complex_double,complex<piranha::double_cf> > ancestor;
      typedef piranha::numerical_container_complex_toolbox<piranha::double_cf> complex_toolbox;
      typedef complex self;
      friend class piranha::numerical_container_complex_toolbox<piranha::double_cf>;
    public:
      typedef piranha::double_cf value_type;
      using ancestor::swap;
      using ancestor::print_plain;
      using ancestor::print_latex;
      using ancestor::is_compatible;
      using ancestor::checkup;
      using ancestor::invert_sign;
      using ancestor::t_eval;
      using ancestor::norm;
      using ancestor::is_zero;
      using ancestor::increase_size;
      using ancestor::append_args;
      using ancestor::prepend_args;
      using ancestor::add;
      using ancestor::subtract;
      using ancestor::mult_by;
      using ancestor::mult_by_self;
      using complex_toolbox::mult_by_self;
      using ancestor::divide_by;
// Start implementation of basic pseries coefficient interface.
//------------
// Basic ctors and dtor.
      explicit complex():ancestor::numerical_container() {}
      explicit complex(const std::string &s):ancestor::numerical_container(s) {}
      explicit complex(const piranha::psymbol &):ancestor::numerical_container() {}
      explicit complex(int n):ancestor::numerical_container(n) {}
      explicit complex(const double &x):ancestor::numerical_container(x) {}
      complex(const complex &c):ancestor::numerical_container(c) {}
      ~complex() {}
// Complex specific contructors.
      explicit complex(int r, int i):complex_toolbox::numerical_container_complex_toolbox(r,i) {}
      explicit complex(const std::complex<int> &c):complex_toolbox::numerical_container_complex_toolbox(c) {}
      explicit complex(const double &r, const double &i):complex_toolbox::numerical_container_complex_toolbox(r,i) {}
      explicit complex(const std::complex<double> &c):complex_toolbox::numerical_container_complex_toolbox(c) {}
      explicit complex(const value_type &r):complex_toolbox::numerical_container_complex_toolbox(r) {}
      explicit complex(const value_type &r, const value_type &i):
        complex_toolbox::numerical_container_complex_toolbox(r,i) {}
// Operators.
      self &operator=(const self &val2)
      {
        return assign_self(val2);
      }
      self &operator=(const value_type &r2)
      {
        return assign_self(r2);
      }
// End implementation of complex basic pseries coefficient interface.
//------------
  };
}
#endif
