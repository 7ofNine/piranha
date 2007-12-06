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

#ifndef PIRANHA_BASIC_PSERIES_COEFFICIENT_CONCEPT_H
#define PIRANHA_BASIC_PSERIES_COEFFICIENT_CONCEPT_H

#include <boost/static_assert.hpp>
#include <string>

#include "../psymbol.h"
#include "../type_traits/eval_type.h"

#define __STATIC_ASSERTION_FAILURE BOOST_STATIC_ASSERT(sizeof(Model) == 0)

namespace piranha
{
/// Concept class for basic coefficients of Poisson series.
/**
 * This concept defines a Poisson series coefficient which enables basic operations on
 * Poisson series.
 */
  template <class Model>
    class basic_pseries_coefficient_concept
  {
    public:
// Ctors and dtor.
/// Default constructor.
// Cannot use static assert here since it will be always called.
      explicit basic_pseries_coefficient_concept() {}
/// Constructor from string.
      explicit basic_pseries_coefficient_concept(const std::string &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from psymbol.
      explicit basic_pseries_coefficient_concept(const psymbol &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from int.
      explicit basic_pseries_coefficient_concept(int)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from double.
      explicit basic_pseries_coefficient_concept(const double &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Copy constructor.
      basic_pseries_coefficient_concept(const Model &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Destructor.
      ~basic_pseries_coefficient_concept() {}
/// Assignment operator.
      Model &operator=(const Model &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Calculate norm.
      double norm(const vector_psym_p &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return double();
      }
/// Is coefficient zero?
      bool is_zero(const vector_psym_p &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return bool();
      }
/// Diagnostic checkup.
// TODO: this should become a templated method.
      bool checkup(const size_t &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return bool();
      }
      bool is_compatible(const size_t &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return bool();
      }
/// Evaluation.
      typename eval_type<basic_pseries_coefficient_concept>::type t_eval(const double &, const vector_psym_p &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return double();
      }
/// Bessel function of the first kind of order n.
// TODO: move into separate concept.
      Model besselJ(int n, const vector_psym_p &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return Model();
      }
/// Print in plain format.
      void print_plain(std::ostream &, const vector_psym_p &) const
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Print in latex format.
      void print_latex(std::ostream &, const vector_psym_p &) const
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Swap.
      Model &swap(Model &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Prepend n arguments.
      void append_args(const size_t &n)
      {
        __STATIC_ASSERTION_FAILURE;
        (void)n;
      }
/// Append n arguments.
      void prepend_args(const size_t &n)
      {
        __STATIC_ASSERTION_FAILURE;
        (void)n;
      }
/// Resize to n arguments.
      void increase_size(const size_t &n)
      {
        __STATIC_ASSERTION_FAILURE;
        (void)n;
      }
// Math ops.
      Model &invert_sign()
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &add(const Model &val2)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &subtract(const Model &val2)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Multiplication with self type.
/**
 * Supports additional templated argument because it may be needed to operate a truncation: the generic template
 * parameter provides a hook from where to fetch information about the truncation method.
 */
      template <class T>
        Model &mult_by_self(const Model &, const T &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &mult_by(int)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &mult_by(const double &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &divide_by(int)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &divide_by(const double &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
  };

/// Concept class for basic complex coefficients of Poisson series.
  template <class RealModel>
    struct complex_basic_pseries_coefficient_concept:
    public basic_pseries_coefficient_concept<std::complex<RealModel> >
  {
    private:
      typedef std::complex<RealModel> Model;
      typedef typename eval_type<Model>::type eval_type;
    public:
// Ctors and dtor.
/// Default constructor.
      explicit complex_basic_pseries_coefficient_concept() {}
/// Constructor from string.
      explicit complex_basic_pseries_coefficient_concept(const std::string &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from psymbol.
      explicit complex_basic_pseries_coefficient_concept(const psymbol &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from int.
      explicit complex_basic_pseries_coefficient_concept(int)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from double.
      explicit complex_basic_pseries_coefficient_concept(const double &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from complex integer.
      explicit complex_basic_pseries_coefficient_concept(const std::complex<int> &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from real and imaginary integer parts.
      explicit complex_basic_pseries_coefficient_concept(int, int)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from complex double.
      explicit complex_basic_pseries_coefficient_concept(const std::complex<double> &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from real and imaginary double parts.
      explicit complex_basic_pseries_coefficient_concept(const double &, const double &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from real type.
      explicit complex_basic_pseries_coefficient_concept(const RealModel &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Constructor from real and imaginary parts.
      explicit complex_basic_pseries_coefficient_concept(const RealModel &, const RealModel &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Copy constructor.
      complex_basic_pseries_coefficient_concept(const Model &)
      {
        __STATIC_ASSERTION_FAILURE;
      }
/// Destructor.
      ~complex_basic_pseries_coefficient_concept() {}
/// Assignment operator.
      Model &operator=(const Model &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Assignment to real type.
      Model &operator=(const RealModel &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Multiply by complex int.
      Model &mult_by(const std::complex<int> &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Multiply by complex double.
      Model &mult_by(const std::complex<double> &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
/// Multiply by real type.
      template <class T>
        Model &mult_by_self(const RealModel &, const T &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
  };
}

#undef __STATIC_ASSERTION_FAILURE

#endif
