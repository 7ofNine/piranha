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

#ifndef PIRANHA_SIMPLE_CONTAINER_H
#define PIRANHA_SIMPLE_CONTAINER_H

#include <iostream>
#include <string>

#include "psymbol_manager.h"
#include "utils.h"            // Lexical converter.

namespace piranha
  {
  /// Simple container class.
  /**
   * This class can be used as a base class for coefficients that consist of a simple entity
   * (float, mpf_class, long double, etc.).
   */
  template <class T>
  class simple_container
    {
    public:
      /// Alias for self.
      typedef simple_container self;
      // Start INTERFACE definition.
      //-------------------------------------------------------
      // Ctors.
      /// Default constructor.
      explicit simple_container():value_(T())
      {}
      /// Constructor from T value.
      explicit simple_container(const T &val):value_(val)
      {}
      /// Copy constructor.
      explicit simple_container(const self &sc):value_(sc.value_)
      {}
      /// Constructor from string.
      explicit simple_container(const std::string &s)
      {
        value_=utils::lexical_converter<T>(s);
      }
      /// Destructor.
      ~simple_container()
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
      /// Print in plain format.
      void print_plain(std::ostream &out_stream, const vector_psym_p &) const
        {
          stream_manager::setup_print(out_stream);
          out_stream << value_;
        }
      /// Print in latex format.
      void print_latex(std::ostream &out_stream, const vector_psym_p &) const
        {
          stream_manager::setup_print(out_stream);
          out_stream << "$" << value_ << "$";
        }
      // Manipulation
      /// Swap values with anothe container.
      void swap(self &dc)
      {
        std::swap(value_,dc.value_);
      }
      /// Add argument.
      void add_arg()
      {}
      /// Resize.
      void resize(const size_t &)
      {}
      // Probing.
      /// Diagnostic checkup.
      bool checkup(const size_t &) const
        {
          return true;
        }
      /// Check whether contained value is larger than size.
      // FIXME: maybe here we should check against 0 size?
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
      bool compatible(const size_t &) const
        {
          return true;
        }
    protected:
      T   value_;
    }
  ;
}

#endif
