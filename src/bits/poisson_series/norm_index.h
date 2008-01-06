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

#ifndef PIRANHA_NORM_INDEX_H
#define PIRANHA_NORM_INDEX_H

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include "../arg_manager.h"

namespace piranha
{
  /// Norm extractor.
  /**
   * Extract norm from term, asserting that arguments have been assigned through piranha::arg_manager.
   */
  template <class Term>
    struct norm_extractor
  {
    typedef double result_type;
    double operator()(const Term &t) const
    {
      p_assert(arg_manager::assigned());
      return t.cf().norm(*arg_manager::cf_args());
    }
  };

  /// Norm-based indices for base_pseries.
  /**
   * This class specifies the following indices to be used in piranha::base_pseries: a hashed index for the
   * identification
   * of terms and a norm-sorted index to discard terms in multiplications. The class is to be used as the I
   * parameter in piranha::base_pseries classes.
   */
  template <class Cf, class Trig, template <class, class> class Term>
    struct norm_index
  {
    typedef boost::multi_index::indexed_by <
      boost::multi_index::ordered_unique <
      boost::multi_index::composite_key <
      Term<Cf, Trig>,
      norm_extractor<Term<Cf, Trig> >,
      boost::multi_index::const_mem_fun < Term<Cf, Trig>, const Trig &,
      &Term<Cf, Trig>::trig >
      >,
      boost::multi_index::composite_key_compare<
      std::greater<double>,
      std::less<Trig>
      >
      >,
      boost::multi_index::hashed_unique<boost::multi_index::identity<Term<Cf, Trig> > >
      > type;
  };
}
#endif
