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

#ifndef PIRANHA_MPQ_CF_H
#define PIRANHA_MPQ_CF_H

#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"

namespace piranha
{
  /// Mpq numerical coefficient.
  /**
   * Arbitrary-size rational coefficient type, to be used as coefficient in piranha::base_series.
   *
   * A set of operators is provided to enable interoperability with basic numerical data types.
   */
  class mpq_cf:public numerical_container<mpq_class,mpq_cf>
  {
      // Alias for the parent class.
      typedef numerical_container<mpq_class,mpq_cf> ancestor;
    public:
      // Start implementation of basic pseries coefficient interface.
      //------------
      // Ctors and dtor.
      /// Empty constructor.
      explicit mpq_cf():ancestor::numerical_container() {}
      /// Constructor from string.
      template <class ArgsTuple>
        explicit mpq_cf(const std::string &s, const ArgsTuple &a):ancestor::numerical_container(s,a)
      {
        // We need to canonicalize when reading from string.
        m_value.canonicalize();
      }
      /// Constructor from integer.
      template <class ArgsTuple>
        explicit mpq_cf(const int &val, const ArgsTuple &a):ancestor::numerical_container(val,a) {}
      /// Constructor from double.
      template <class ArgsTuple>
        explicit mpq_cf(const double &val, const ArgsTuple &a):ancestor::numerical_container(val,a) {}
      /// Constructor from psym.
      template <class ArgsTuple>
        explicit mpq_cf(const psym_p &p, const int &n, const ArgsTuple &a):ancestor::numerical_container(p,n,a) {}
      // Override norm and evaluation.
      template <class ArgsTuple>
        double norm(const ArgsTuple &) const
      {
        return std::abs(g_value().get_d());
      }
      // Override division to catch divide by zero.
      template <class ArgsTuple>
        mpq_cf &divide_by(const int &n, const ArgsTuple &a) throw(division_by_zero)
      {
        if (n == 0)
        {
          throw division_by_zero();
        }
        return ancestor::divide_by(n,a);
      }
      template <class ArgsTuple>
        mpq_cf &divide_by(const double &x, const ArgsTuple &a) throw(division_by_zero)
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
        return g_value().get_d();
      }
      int get_int() const throw (unsuitable)
      {
        const char *msg = "Cannot convert rational coefficient to integer.";
        int retval = ancestor::m_value.get_num().get_si();
        if (ancestor::m_value.get_den() != 1)
        {
          throw (unsuitable(msg));
        }
        return retval;
      }
  };
}

#endif
