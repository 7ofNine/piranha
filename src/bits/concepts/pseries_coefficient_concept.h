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

#ifndef PIRANHA_PSERIES_COEFFICIENT_CONCEPT_H
#define PIRANHA_PSERIES_COEFFICIENT_CONCEPT_H

#include <boost/static_assert.hpp>
#include <complex>
#include <iostream>

#include "../p_assert.h"
#include "../psymbol.h"
#include "../type_traits/eval_type.h"

#define __STATIC_ASSERTION_FAILURE BOOST_STATIC_ASSERT(sizeof(Model) == 0)

namespace piranha
{
// FIXME: deal with operator=() !!.
/// Concept class for coefficients of Poisson series.
  template <class Model>
    class pseries_coefficient_concept
  {
    public:
// Ctors and dtor.
/// Default constructor.
      explicit pseries_coefficient_concept()
      {}
/// Constructor from string.
      explicit pseries_coefficient_concept(const std::string &)
      {
        hard_assert(false);
      }
/// Constructor from double.
      explicit pseries_coefficient_concept(const double &)
      {
        hard_assert(false);
      }
/// Constructor from integer.
      explicit pseries_coefficient_concept(int)
      {
        hard_assert(false);
      }
/// Constructor from symbol.
      explicit pseries_coefficient_concept(const psymbol &)
      {
        hard_assert(false);
      }
/// Copy constructor.
      pseries_coefficient_concept(const Model &)
      {
        hard_assert(false);
      }
/// Destructor.
      ~pseries_coefficient_concept()
      {}
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
      bool compatible(const size_t &) const
      {
        __STATIC_ASSERTION_FAILURE;
        return bool();
      }
/// Evaluation.
      typename eval_type<pseries_coefficient_concept>::type t_eval(const double &, const vector_psym_p &) const
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
// Math ops.
      Model &invert_sign()
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &add_self(const Model &val2)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &subtract_self(const Model &val2)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
// TODO: use proper doxygen lists in documentation below.
/// Multiplication with self type.
/**
 * Supports additional templated arguments because:
 *
 * - we want interoperability complex <--> real;
 *
 * - it may be needed to operate a truncation: the second generic template
 * parameter provides a hook from where to fetch information about the truncation method.
 */
      template <class T, class U>
        Model &mult_by_self(const T &, const U &)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &mult_by_int(int n)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &mult_by_double(const double &x)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      template <class T>
        Model &mult_by_generic(const T &x)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      Model &divide_by_int(int n)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
      template <class T>
        Model &divide_by_generic(const T &x)
      {
        __STATIC_ASSERTION_FAILURE;
        return *static_cast<Model *>(this);
      }
  };
}

namespace std
{
/// Concept class for complex coefficients of Poisson series.
/**
 * Partial template specialization of standard complex class.
 */
  template <>
    template <class RealModel>
    class complex<piranha::pseries_coefficient_concept<RealModel> >:
    public piranha::pseries_coefficient_concept<complex<RealModel> >
  {
      typedef complex<RealModel> Model;
      typedef typename piranha::eval_type<complex<Model> >::type eval_type;
    public:
// Ctors and dtor.
/// Default constructor.
      explicit complex()
      {}
/// Constructor from string.
      explicit complex(const std::string &)
      {
        hard_assert(false);
      }
/// Constructor from double.
      explicit complex(const double &)
      {
        hard_assert(false);
      }
/// Constructor from integer.
      explicit complex(int)
      {
        hard_assert(false);
      }
/// Constructor from symbol.
      explicit complex(const piranha::psymbol &)
      {
        hard_assert(false);
      }
/// Copy constructor.
      complex(const Model &)
      {
        hard_assert(false);
      }
/// Destructor.
      ~complex()
      {}
  };
}

#undef __STATIC_ASSERTION_FAILURE

#endif
