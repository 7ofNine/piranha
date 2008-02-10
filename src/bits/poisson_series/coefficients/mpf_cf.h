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

#ifndef PIRANHA_MPF_CF_H
#define PIRANHA_MPF_CF_H

#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../../common_typedefs.h"
#include "../../math.h"                           // besselJ.
#include "../../psymbol.h"
#include "../base_classes/numerical_container.h"
#include "../concepts/basic_pseries_coefficient_concept.h"

namespace piranha
{
  /// Mpz numerical coefficient.
  /**
   * This class can be used as coefficient in Poisson series. It encapsulate a GMP mpf_class
   * value and provides access to it.
   *
   * A set of operators is provided to enable interoperability with basic numerical data types.
   */
  class mpf_cf:
  public concepts::basic_pseries_coefficient_concept<mpf_cf>,
    public numerical_container<mpf_class,mpf_cf>
  {
    /// Alias for self.
    typedef mpf_cf self;
    /// Alias for the parent class.
    typedef numerical_container<mpf_class,mpf_cf> ancestor;
    /// Alias for concept.
    typedef concepts::basic_pseries_coefficient_concept<mpf_cf> concept;
    public:
      // This is necessary because some moethods are present in concept _and_ in numerical container.
      // We avoid the dreaded diamond problem by explicitly stating which functions to use.
      using ancestor::swap;
      using ancestor::print_plain;
      using ancestor::print_latex;
      using ancestor::checkup;
      using ancestor::invert_sign;
      using ancestor::is_ignorable;
      using ancestor::is_insertable;
      using ancestor::needs_padding;
      using ancestor::pad_right;
      using ancestor::apply_layout;
      using ancestor::add;
      using ancestor::subtract;
      using ancestor::mult_by;
      using ancestor::mult_by_self;
      using ancestor::divide_by;
      // Start implementation of basic pseries coefficient interface.
      //------------
      // Ctors and dtor.
      /// Empty constructor.
      explicit mpf_cf():ancestor::numerical_container() {}
      /// Constructor from string.
      explicit mpf_cf(const std::string &s):ancestor::numerical_container(s) {}
      /// Constructor from symbol.
      explicit mpf_cf(const psymbol &):ancestor::numerical_container()
      {
        std::cout << "WARNING: building numerical coefficient from psymbol." << std::endl;
      }
      /// Constructor from integer.
      explicit mpf_cf(int val):ancestor::numerical_container(val) {}
      /// Constructor from double.
      explicit mpf_cf(const double &val):ancestor::numerical_container(val) {}
      /// Copy constructor.
      mpf_cf(const self &dc):concept::basic_pseries_coefficient_concept(dc),ancestor::numerical_container(dc) {}
      /// Destructor.
      ~mpf_cf() {}
      // Needed operators.
      self &operator=(const self &val2)
      {
        return assign_self(val2);
      }
      // Override norm and evaluation.
      template <class ArgsTuple>
        double norm(const ArgsTuple &) const
      {
        return std::abs(g_value().get_d());
      }
      double t_eval(const double &, const vector_psym_p &) const
      {
        return g_value().get_d();
      }
      static const size_t max_size = ancestor::max_size;
      // End implementation of basic pseries coefficient interface.
      //------------
      // Start implementation of trigonometric pseries coefficient interface.
      // Used in:
      // - trigonometric toolbox,
      //------------
      // TODO: Move into own toolbox and concept.
      self besselJ(int n, const vector_psym_p &) const
      {
        self retval;
        retval.s_value()=math::besselJ(n,g_value().get_d());
        return retval;
      }
      // End implementation of trigonometric pseries coefficient interface.
      //------------
      // Start implementation of power-enabled pseries coefficient interface.
      // TODO: Move into own toolbox and concept.
      self pow(const double &y) const
      {
        self retval;
        retval.s_value()=std::pow(g_value().get_d(),y);
        return retval;
      }
  };
}


namespace std
{
  template <>
    struct complex<piranha::mpf_cf>:
  public piranha::concepts::complex_basic_pseries_coefficient_concept<piranha::mpf_cf>,
    public piranha::numerical_container<complex<mpf_class>,complex<piranha::mpf_cf> >,
    public piranha::numerical_container_complex_toolbox<piranha::mpf_cf>
  {
    private:
      typedef piranha::numerical_container<complex<mpf_class>,complex<piranha::mpf_cf> > ancestor;
      typedef piranha::numerical_container_complex_toolbox<piranha::mpf_cf> complex_toolbox;
      typedef complex self;
      typedef piranha::concepts::complex_basic_pseries_coefficient_concept<piranha::mpf_cf> concept;
      friend class piranha::numerical_container_complex_toolbox<piranha::mpf_cf>;
    public:
      typedef piranha::mpf_cf value_type;
      using ancestor::swap;
      using ancestor::print_plain;
      using ancestor::print_latex;
      using ancestor::checkup;
      using ancestor::invert_sign;
      using ancestor::is_ignorable;
      using ancestor::is_insertable;
      using ancestor::needs_padding;
      using ancestor::pad_right;
      using ancestor::apply_layout;
      using ancestor::add;
      using ancestor::subtract;
      using ancestor::mult_by;
      using ancestor::mult_by_self;
      using ancestor::divide_by;
      using complex_toolbox::mult_by;
      using complex_toolbox::divide_by;
      using complex_toolbox::mult_by_self;
      using complex_toolbox::real;
      using complex_toolbox::imag;
      using complex_toolbox::set_real;
      using complex_toolbox::set_imag;
      // Start implementation of basic pseries coefficient interface.
      //------------
      // Basic ctors and dtor.
      explicit complex():ancestor::numerical_container() {}
      explicit complex(const std::string &s):ancestor::numerical_container(s) {}
      explicit complex(const piranha::psymbol &):ancestor::numerical_container() {}
      explicit complex(int n):ancestor::numerical_container(n) {}
      explicit complex(const double &x):ancestor::numerical_container(x) {}
      complex(const complex &c):concept::complex_basic_pseries_coefficient_concept(c),
        ancestor::numerical_container(c),complex_toolbox::numerical_container_complex_toolbox(c) {}
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
      // Probing.
      double norm(const piranha::vector_psym_p &) const
      {
        // NOTICE: the success of this probably depends upon std::complex implementation...
        return std::abs(g_value()).get_d();
      }
      // Override evaluation.
      std::complex<double> t_eval(const double &, const piranha::vector_psym_p &) const
      {
        return std::complex<double>(g_value().real().get_d(),g_value().imag().get_d());
      }
      static const size_t max_size = ancestor::max_size;
      // End implementation of complex basic pseries coefficient interface.
      //------------
  };
}
#endif
