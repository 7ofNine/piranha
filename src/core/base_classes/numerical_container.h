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

#ifndef PIRANHA_NUMERICAL_CONTAINER_H
#define PIRANHA_NUMERICAL_CONTAINER_H

#include <iostream>
#include <string>

#include "../arg_manager.h"
#include "../psym.h"
#include "../utils.h" // Lexical converter.
#include "../type_traits.h"

// Convenience macros.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  /// Numerical container class.
  /**
   * This class can be used as a base class for coefficients that consist of a
   * numerical entity (double, GMP classes, etc.).
   */
  template <class T, class Derived>
    class numerical_container
  {
      /// Alias for evaluation type.
      typedef typename eval_type<Derived>::type eval_type;
      /// Alias for self.
      typedef numerical_container self;
    public:
      // Start implementation of basic pseries coefficient interface.
      //------------
      // Ctors.
      explicit numerical_container():m_value(0) {}
      template <class ArgsTuple>
        explicit numerical_container(const std::string &s, const ArgsTuple &):m_value(utils::lexical_converter<T>(s))
      {}
      template <class ArgsTuple>
        explicit numerical_container(const int &n, const ArgsTuple &):m_value(n) {}
      template <class ArgsTuple>
        explicit numerical_container(const double &x, const ArgsTuple &):m_value(x) {}
      /// Ctor from psym.
      /**
       * Sets m_value to one.
       */
      template <class ArgsTuple>
        explicit numerical_container(const psym_p &, const int &, const ArgsTuple &):m_value(1) {}
      // I/O.
      template <class ArgsTuple>
        void print_plain(std::ostream &out_stream, const ArgsTuple &) const
      {
        stream_manager::setup_print(out_stream);
        out_stream << m_value;
      }
      template <class ArgsTuple>
      void print_latex(std::ostream &out_stream, const ArgsTuple &) const
      {
// TODO: rework this.
//         stream_manager::setup_print(out_stream);
//         out_stream << "$" << m_value << "$";
      }
      // Manipulation
      Derived &swap(Derived &dc)
      {
        std::swap(m_value,dc.m_value);
        return *derived_cast;
      }
      template <class ArgsTuple>
        void pad_right(const ArgsTuple &)
      {}
      template <class ArgsTuple, class Layout>
        void apply_layout(const ArgsTuple &, const Layout &)
      {}
      // Probing.
      template <class ArgsTuple>
        bool checkup(const ArgsTuple &) const
      {
        return true;
      }
      template <class ArgsTuple>
        double norm(const ArgsTuple &) const
      {
        return std::abs(m_value);;
      }
      static const size_t max_size = 0;
      // If value is less than numericalzero  in absolute value it is considered
      // to be zero.
      template <class ArgsTuple>
        bool is_ignorable(const ArgsTuple &a) const
      {
        return (static_cast<Derived const *>(this)->norm(a) < settings::numerical_zero());
      }
      template <class ArgsTuple>
        bool is_insertable(const ArgsTuple &) const
      {
        return true;
      }
      template <class ArgsTuple>
        bool needs_padding(const ArgsTuple &) const
      {
        return false;
      }
      template <class ArgsTuple>
        const eval_type &eval(const double &, const ArgsTuple &) const
      {
        return m_value;
      }
      // Maths.
      template <class ArgsTuple>
        void invert_sign(const ArgsTuple &)
      {
        m_value*=-1;
      }
      template <class ArgsTuple>
        Derived &add(const Derived &val2, const ArgsTuple &)
      {
        return add_generic(val2.g_value());
      }
      template <class ArgsTuple>
        Derived &subtract(const Derived &val2, const ArgsTuple &)
      {
        return subtract_generic(val2.g_value());
      }
      template <class ArgsTuple>
        Derived &mult_by(const int &n, const ArgsTuple &)
      {
        return mult_by_generic(n);
      }
      template <class ArgsTuple>
        Derived &mult_by(const double &x, const ArgsTuple &)
      {
        return mult_by_generic(x);
      }
      template <class ArgsTuple>
        Derived &mult_by(const self &x, const ArgsTuple &)
      {
        return mult_by_generic(x.g_value());
      }
      template <class ArgsTuple>
        Derived &divide_by(const int &n, const ArgsTuple &)
      {
        return divide_by_generic(n);
      }
      template <class ArgsTuple>
        Derived &divide_by(const double &x, const ArgsTuple &)
      {
        return divide_by_generic(x);
      }
      // Multiply and add.
      template <class ArgsTuple>
        void addmul(const Derived &x1, const Derived &x2, const ArgsTuple &)
      {
        m_value += x1.m_value * x2.m_value;
      }
      // End implementation of basic pseries coefficient interface.
      //------------
      // TODO: move into own toolbox.
      void partial(const size_t &, Derived &retval) const
      {
        retval=Derived(0);
      }
      /// Get value.
      const T &g_value() const
      {
        return m_value;
      }
      /// Set value.
      T &s_value()
      {
        return m_value;
      }
    protected:
      template <class U>
        Derived &assign_self(const U &x)
      {
        m_value=x.g_value();
        return *derived_cast;
      }
      template <class U>
        Derived &add_generic(const U &x)
      {
        m_value+=x;
        return *derived_cast;
      }
      template <class U>
        Derived &subtract_generic(const U &x)
      {
        m_value-=x;
        return *derived_cast;
      }
      template <class U>
        Derived &mult_by_generic(const U &x)
      {
        m_value*=x;
        return *derived_cast;
      }
      template <class U>
        Derived &divide_by_generic(const U &x)
      {
        m_value/=x;
        return *derived_cast;
      }
    protected:
      // Data member.
      T m_value;
  };

  /// Toolbox for complex-specific methods of piranha::numerical_container.
  template <class realDerived>
    class numerical_container_complex_toolbox
  {
      typedef std::complex<realDerived> Derived;
      typedef realDerived value_type;
    public:
      // TODO: drop usageof g_value() - s_value() to improve performance?
      // Start implementation of complex basic pseries coefficient interface.
      //------------
      // Ctors.
      numerical_container_complex_toolbox() {}
      explicit numerical_container_complex_toolbox(int r, int i)
      {
        derived_cast->s_value().real()=r;
        derived_cast->s_value().imag()=i;
      }
      explicit numerical_container_complex_toolbox(const std::complex<int> &c)
      {
        derived_cast->s_value()=c;
      }
      explicit numerical_container_complex_toolbox(const double &r, const double &i)
      {
        derived_cast->s_value().real()=r;
        derived_cast->s_value().imag()=i;
      }
      explicit numerical_container_complex_toolbox(const std::complex<double> &c)
      {
        derived_cast->s_value()=c;
      }
      explicit numerical_container_complex_toolbox(const value_type &r)
      {
        derived_cast->s_value().real()=r.g_value();
      }
      explicit numerical_container_complex_toolbox(const value_type &r, const value_type &i)
      {
        derived_cast->s_value().real()=r.g_value();
        derived_cast->s_value().imag()=i.g_value();
      }
      // Getters and setters.
      value_type real() const
      {
        value_type retval;
        retval.s_value()=derived_cast->g_value().real();
        return retval;
      }
      value_type imag() const
      {
        value_type retval;
        retval.s_value()=derived_cast->g_value().imag();
        return retval;
      }
      void set_real(const value_type &r)
      {
        derived_cast->s_value()=r.g_value();
      }
      void set_imag(const value_type &i)
      {
        derived_cast->s_value().real()=0;
        derived_cast->s_value().imag()=i.g_value();
      }
      // Maths.
      template <class DerivedPs>
        Derived &mult_by_self(const value_type &x, const DerivedPs &)
      {
        derived_cast->mult_by_generic(x.g_value());
      }
      Derived &mult_by(const std::complex<int> &c)
      {
        derived_cast->mult_by_generic(c);
      }
      Derived &mult_by(const std::complex<double> &c)
      {
        derived_cast->mult_by_generic(c);
      }
      Derived &divide_by(const std::complex<int> &c)
      {
        derived_cast->divide_by_generic(c);
      }
      Derived &divide_by(const std::complex<double> &c)
      {
        derived_cast->divide_by_generic(c);
      }
      // End implementation of complex basic pseries coefficient interface.
      //------------
  };
}

namespace std
{
  // Overloads for I/O operators.
  template <class T, class Derived>
    inline istream &operator>>(istream &is, piranha::numerical_container<T,Derived> &nc)
  {
    string tmp;
    getline(is,tmp);
    nc.s_value()=piranha::utils::lexical_converter<T>(tmp);
    return is;
  }

  template <class T, class Derived>
    inline ostream &operator<<(ostream &os, const piranha::numerical_container<T,Derived> &nc)
  {
    os << nc.g_value();
    return os;
  }
}

#undef derived_const_cast
#undef derived_cast

#endif