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

#ifndef PIRANHA_NUMERICAL_CONTAINER_H
#define PIRANHA_NUMERICAL_CONTAINER_H

#include <iostream>
#include <string>

#include "../bits/psymbol.h"
#include "../bits/utils.h"                                // Lexical converter.
#include "../bits/type_traits/eval_type.h"

namespace piranha
{
/// Simple container class.
/**
 * This class can be used as a base class for coefficients that consist of a
 * numerical entity (double, GMP classes, etc.).
 */
  template <class T, class Derived>
    class numerical_container
  {
      typedef typename eval_type<Derived>::type eval_type;
    public:
/// Alias for self.
      typedef numerical_container self;
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default constructor.
      explicit numerical_container():value_(0)
        {}
/// Constructor from T value.
      explicit numerical_container(const T &val):value_(val)
        {}
/// Copy constructor.
      explicit numerical_container(const self &sc):value_(sc.value_)
        {}
/// Constructor from string.
      explicit numerical_container(const std::string &s)
      {
        value_=utils::lexical_converter<T>(s);
      }
/// Destructor.
      ~numerical_container()
        {}
// Getters.
/// Get value.
      const T &value() const
      {
        return value_;
      }
/// Get actual width.
      size_t actual_width() const
      {
        return 0;
      }
// Setters
/// Set value.
      T &value()
      {
        return value_;
      }
// I/O.
      void print_plain(std::ostream &out_stream, const vector_psym_p &) const
      {
        stream_manager::setup_print(out_stream);
        out_stream << value_;
      }
      void print_latex(std::ostream &out_stream, const vector_psym_p &) const
      {
        stream_manager::setup_print(out_stream);
        out_stream << "$" << value_ << "$";
      }
// Manipulation
      Derived &swap(Derived &dc)
      {
        std::swap(value_,dc.value_);
        return *static_cast<Derived *>(this);
      }
/// Prepend arguments.
// TODO: place asserts here, to check we never want to resize to > 0.
// TODO: how does this interact with appending arguments from series?
// The problem here is how to handle resize request. Maybe coefficient
// and trigs should have a trait that tells whether they are resizable or not?
      void append_args(const size_t &)
        {}
/// Append arguments.
      void prepend_args(const size_t &)
        {}
/// Resize.
      void increase_size(const size_t &)
        {}
// Probing.
/// Diagnostic checkup.
      bool checkup(const size_t &) const
      {
        return true;
      }
/// Check whether contained value is larger than size.
// TODO: maybe here we should check against 0 size?
      bool larger(const size_t &) const
      {
        return false;
      }
/// Check whether contained value is smaller than size.
      bool smaller(const size_t &) const
      {
        return false;
      }
/// Check whether contained value is size compatible.
      bool is_compatible(const size_t &) const
      {
        return true;
      }
/// Evaluation.
/**
 * Evaluation for this class always returns the same value.
 */
      const eval_type &t_eval(const double &, const vector_psym_p &) const
      {
        return value_;
      }
/// Partial derivative.
/**
 * Always returns 0, since this is a purely numerical quantity.
 * @param[out] retval, Derived return value.
 */
      void partial(const size_t &, Derived &retval) const
      {
        retval=Derived(0);
      }
// Maths.
      Derived &invert_sign()
      {
        value_*=-1;
        return *static_cast<Derived *>(this);
      }
      Derived &add_self(const Derived &val2)
      {
        value_+=val2.value();
        return *static_cast<Derived *>(this);
      }
      Derived &subtract_self(const Derived &val2)
      {
        value_-=val2.value();
        return *static_cast<Derived *>(this);
      }
      self &mult_by_int(int n)
      {
        return mult_by_generic(n);
      }
      self &mult_by_double(const double &x)
      {
        return mult_by_generic(x);
      }
      template <class U>
        self &mult_by_generic(const U &x)
      {
        value_*=x;
        return *static_cast<Derived *>(this);
      }
      template <class U, class DerivedPs>
        self &mult_by_self(const U &x, const DerivedPs &)
      {
        value_*=x.value();
        return *static_cast<Derived *>(this);
      }
      self &divide_by_int(int n)
      {
        return divide_by_generic(n);
      }
      template <class U>
        self &divide_by_generic(const U &x)
      {
        value_/=x;
        return *static_cast<Derived *>(this);
      }
    protected:
      T   value_;
  };
}
#endif
