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

#ifndef PIRANHA_BASIC_PSERIES_COEFFICIENT_CONCEPT_H
#define PIRANHA_BASIC_PSERIES_COEFFICIENT_CONCEPT_H

#include <boost/static_assert.hpp>
#include <string>

#include "../../common_typedefs.h"                // For layout_type.
#include "../../psymbol.h"
#include "../../type_traits/eval_type.h"

#define __STATIC_ASSERTION_FAILURE BOOST_STATIC_ASSERT(sizeof(Model) == 0)

namespace piranha
{
  /// Namespace for concept classes.
  /**
   * The classes in this namespace are not to be used, they serve only for documentation purposes.
   */
  namespace concepts
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
        basic_pseries_coefficient_concept(const Model &) {}
        /// Destructor.
        ~basic_pseries_coefficient_concept() {}
        /// Assignment operator.
        Model &operator=(const Model &)
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Calculate norm.
        // TODO: how is this defined for, e.g., polynomials? Abs of evaluation at t=0?
        double norm(const vector_psym_p &) const
        {
          __STATIC_ASSERTION_FAILURE;
          return double();
        }
        /// Is coefficient ignorable?
        template <class Series>
          bool is_ignorable(const Series &) const
        {
          __STATIC_ASSERTION_FAILURE;
          return bool();
        }
        /// Diagnostic checkup.
        // TODO: document why, here and above, we template this way.
        template <class Series>
          bool checkup(const Series &) const
        {
          __STATIC_ASSERTION_FAILURE;
          return bool();
        }
        /// Is coefficient insertable?
        /**
         * Checks whether a term containing the coefficient can be inserted in a series with coefficient width equal
         * to w.
         * @param[in] w size_t coefficient width against which insertability is checked.
         */
        bool is_insertable(const size_t &w) const
        {
          __STATIC_ASSERTION_FAILURE;
          (void)w;
          return bool();
        }
        /// Check whether coefficient needs padding to be inserted in series with coefficient width w.
        bool needs_padding(const size_t &w) const
        {
          __STATIC_ASSERTION_FAILURE;
          (void)w;
          return bool();
        }
        /// Evaluation.
        typename eval_type<basic_pseries_coefficient_concept>::type t_eval(const double &, const vector_psym_p &) const
        {
          __STATIC_ASSERTION_FAILURE;
          return double();
        }
        /// Max size.
        static const size_t max_size = 0;
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
        /// Apply layout.
        void apply_layout(const layout_type &)
        {
          __STATIC_ASSERTION_FAILURE;
        }
        /// Append n arguments and initialise their values to 0.
        void pad_right(const size_t &n)
        {
          __STATIC_ASSERTION_FAILURE;
          (void)n;
        }
        /// Invert sign.
        Model &invert_sign()
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Add self.
        Model &add(const Model &val2)
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Subtract self.
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
        /// Multiply by int.
        Model &mult_by(int)
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Multiply by double.
        Model &mult_by(const double &)
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Divide by int.
        Model &divide_by(int)
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Divide by double.
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
        complex_basic_pseries_coefficient_concept(const Model &) {}
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
        /// Divide by complex int.
        Model &divide_by(const std::complex<int> &)
        {
          __STATIC_ASSERTION_FAILURE;
          return *static_cast<Model *>(this);
        }
        /// Divide by complex double.
        Model &divide_by(const std::complex<double> &)
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
        /// Get copy of real part.
        RealModel real() const
        {
          __STATIC_ASSERTION_FAILURE;
          return RealModel();
        }
        /// Get copy of imaginary part.
        RealModel imag() const
        {
          __STATIC_ASSERTION_FAILURE;
          return RealModel();
        }
        /// Set real part to r.
        void set_real(const RealModel &r)
        {
          __STATIC_ASSERTION_FAILURE;
          (void)r;
        }
        /// Set imaginary part to i.
        void set_imag(const RealModel &i)
        {
          __STATIC_ASSERTION_FAILURE;
          (void)i;
        }
    };
  }
}


#undef __STATIC_ASSERTION_FAILURE
#endif
