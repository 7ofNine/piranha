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

#ifndef PIRANHA_POLYNOMIAL_COMMONS_H
#define	PIRANHA_POLYNOMIAL_COMMONS_H

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include "../arg_manager.h"
#include "../ntuple.h"

namespace piranha
{
  /// Degree extractor.
  /**
   * Extract degree from monomial.
   */
  template <class Term>
    struct degree_extractor
  {
    typedef int result_type;
    int operator()(const Term &m) const
    {
      p_assert((arg_manager<Term>::assigned()));
      return m.m_key.get_degree();
    }
  };

  /// Degree-based index for polynomials.
  /**
   * This class specifies the following indices to be used in polynomials: a degree-sorted index
   * and an exponent-hashed index.
   */
  template <class Monomial>
    struct degree_index
  {
    typedef boost::multi_index::indexed_by <
      boost::multi_index::ordered_non_unique <
      degree_extractor<Monomial>
      >,
      boost::multi_index::hashed_unique<boost::multi_index::identity<Monomial> >
      > type;
  };
}

#endif

