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
#include "../p_assert.h"
#include "../psymbol.h"

#define __MODEL_SIZE_STATIC_ASSERTION BOOST_STATIC_ASSERT(sizeof(Model) == 0)

namespace piranha
{
/// Concept class for coefficients of Poisson series.
  template <class Model>
    class pseries_coefficient_concept
  {
    public:
// Ctors and dtor.
/// Default constructor.
      explicit pseries_coefficient_concept()
      {}
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
/// Constructor from string.
      explicit pseries_coefficient_concept(const std::string &)
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
        __MODEL_SIZE_STATIC_ASSERTION;
        return double();
      }
/// Is coefficient zero?
      bool is_zero(const vector_psym_p &) const
      {
        __MODEL_SIZE_STATIC_ASSERTION;
        return bool();
      }
/// Bessel function of the first kind of order n.
// TODO: move into separate concept.
      Model besselJ(int n, const vector_psym_p &) const
      {
        __MODEL_SIZE_STATIC_ASSERTION;
        return Model();
      }
  };
}

#undef __MODEL_SIZE_STATIC_ASSERTION

#endif
