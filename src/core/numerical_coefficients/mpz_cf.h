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

#ifndef PIRANHA_MPZ_CF_H
#define PIRANHA_MPZ_CF_H

#include <cmath>
#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../settings.h"

namespace piranha
{
  /// Mpz numerical coefficient.
  /**
   * Arbitrary-size integer coefficient type, to be used as coefficient in piranha::base_series.
   *
   * A set of operators is provided to enable interoperability with basic numerical data types.
   */
  class mpz_cf:public numerical_container<mpz_class,mpz_cf>
  {
      // Alias for the parent class.
      typedef numerical_container<mpz_class,mpz_cf> ancestor;
    public:
      // Ctors and dtor.
      /// Empty constructor.
      explicit mpz_cf():ancestor::numerical_container() {}
      /// Constructor from string.
      template <class ArgsTuple>
        explicit mpz_cf(const std::string &s, const ArgsTuple &a):ancestor::numerical_container(s,a)
      {}
      /// Constructor from integer.
      template <class ArgsTuple>
        explicit mpz_cf(const max_fast_int &val, const ArgsTuple &a):ancestor::numerical_container(val,a) {}
      /// Constructor from double.
      template <class ArgsTuple>
        explicit mpz_cf(const double &val, const ArgsTuple &a):ancestor::numerical_container(val,a) {}
      /// Constructor from psym.
      template <class ArgsTuple>
        explicit mpz_cf(const psym_p &p, const int &n, const ArgsTuple &a):ancestor::numerical_container(p,n,a) {}
      // Override norm and evaluation.
      template <class ArgsTuple>
        double norm(const ArgsTuple &) const
      {
        return std::abs(g_value().get_d());
      }
      // Override division to catch divide by zero.
      template <class ArgsTuple>
        mpz_cf &divide_by(const max_fast_int &n, const ArgsTuple &a) throw(division_by_zero)
      {
        if (n == 0)
        {
          throw division_by_zero();
        }
        return ancestor::divide_by(n,a);
      }
      template <class ArgsTuple>
        mpz_cf &divide_by(const double &x, const ArgsTuple &a) throw(division_by_zero)
      {
        if (x == 0)
        {
          throw division_by_zero();
        }
        return ancestor::divide_by(x,a);
      }
      // Override this, hence avoiding to calculate norm.
      template <class ArgsTuple>
        bool is_ignorable(const ArgsTuple &) const
      {
        return (m_value == 0);
      }
      template <class ArgsTuple>
        double eval(const double &, const ArgsTuple &) const
      {
        return m_value.get_d();
      }
      // Multiply and add.
      template <class ArgsTuple>
        void addmul(const mpz_cf &x1, const mpz_cf &x2, const ArgsTuple &)
      {
        mpz_addmul(m_value.get_mpz_t(),x1.m_value.get_mpz_t(),x2.m_value.get_mpz_t());
      }
      template <class ArgsTuple>
        mpz_cf pow(const double &y, const ArgsTuple &) const
      {
        mpz_cf retval;
        // If value = 1, then any power is ok, just return 1.
        if (m_value == 1)
        {
          retval.m_value = 1;
          return retval;
        }
        const max_fast_int pow_n((max_fast_int)nearbyint(y));
        if (std::abs(pow_n - y) > settings::numerical_zero())
        {
          throw (unsuitable("Cannot raise integer coefficient different from unity to real power."));
        }
        if (pow_n < 0)
        {
          throw (unsuitable("Cannot raise integer coefficient different from unity to negative integer power."));
        }
        mpz_pow_ui(retval.m_value.get_mpz_t(),m_value.get_mpz_t(),(size_t)pow_n);
        return retval;
      }
  };
}


namespace std
{
//   template <>
//     struct complex<piranha::mpf_cf>:
//     public piranha::numerical_container<complex<mpf_class>,complex<piranha::mpf_cf> >,
//     public piranha::numerical_container_complex_toolbox<piranha::mpf_cf>
//   {
//     private:
//       typedef piranha::numerical_container<complex<mpf_class>,complex<piranha::mpf_cf> > ancestor;
//       typedef piranha::numerical_container_complex_toolbox<piranha::mpf_cf> complex_toolbox;
//       typedef complex self;
//       friend class piranha::numerical_container_complex_toolbox<piranha::mpf_cf>;
//     public:
//       typedef piranha::mpf_cf value_type;
// //       using ancestor::swap;
// //       using ancestor::print_plain;
// //       using ancestor::print_latex;
// //       using ancestor::checkup;
// //       using ancestor::invert_sign;
// //       using ancestor::is_ignorable;
// //       using ancestor::is_insertable;
// //       using ancestor::needs_padding;
// //       using ancestor::pad_right;
// //       using ancestor::apply_layout;
// //       using ancestor::add;
// //       using ancestor::subtract;
// //       using ancestor::mult_by;
// //       using ancestor::mult_by_self;
// //       using ancestor::divide_by;
// //       using complex_toolbox::mult_by;
// //       using complex_toolbox::divide_by;
// //       using complex_toolbox::mult_by_self;
// //       using complex_toolbox::real;
// //       using complex_toolbox::imag;
// //       using complex_toolbox::set_real;
// //       using complex_toolbox::set_imag;
//       // Start implementation of basic pseries coefficient interface.
//       //------------
//       // Basic ctors and dtor.
//       explicit complex():ancestor::numerical_container() {}
//       template <class ArgsTuple>
//         explicit complex(const std::string &s, const ArgsTuple &a):ancestor::numerical_container(s,a) {}
//       explicit complex(int n):ancestor::numerical_container(n) {}
//       explicit complex(const double &x):ancestor::numerical_container(x) {}
//       complex(const complex &c):ancestor::numerical_container(c),complex_toolbox::numerical_container_complex_toolbox(c) {}
//       // Complex specific contructors.
//       explicit complex(int r, int i):complex_toolbox::numerical_container_complex_toolbox(r,i) {}
//       explicit complex(const std::complex<int> &c):complex_toolbox::numerical_container_complex_toolbox(c) {}
//       explicit complex(const double &r, const double &i):complex_toolbox::numerical_container_complex_toolbox(r,i) {}
//       explicit complex(const std::complex<double> &c):complex_toolbox::numerical_container_complex_toolbox(c) {}
//       explicit complex(const value_type &r):complex_toolbox::numerical_container_complex_toolbox(r) {}
//       explicit complex(const value_type &r, const value_type &i):
//       complex_toolbox::numerical_container_complex_toolbox(r,i) {}
//       // Operators.
//       self &operator=(const self &val2)
//       {
//         return assign_self(val2);
//       }
//       self &operator=(const value_type &r2)
//       {
//         return assign_self(r2);
//       }
//       // Probing.
//       template <class ArgsTuple>
//         double norm(const ArgsTuple &) const
//       {
//         // NOTICE: the success of this probably depends upon std::complex implementation...
//         return std::abs(g_value()).get_d();
//       }
//       // Override evaluation.
//       template <class ArgsTuple>
//         std::complex<double> t_eval(const double &, const ArgsTuple &) const
//       {
//         return std::complex<double>(g_value().real().get_d(),g_value().imag().get_d());
//       }
//       // End implementation of complex basic pseries coefficient interface.
//       //------------
//   };
}
#endif
