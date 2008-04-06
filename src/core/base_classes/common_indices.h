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

#ifndef PIRANHA_COMMON_INDICES_H
#define PIRANHA_COMMON_INDICES_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include "../arg_manager.h"

namespace piranha
{
  // Key extractor.
  template <class Term>
    struct key_extractor
  {
    typedef const typename Term::key_type &result_type;
    const typename Term::key_type &operator()(const Term &t) const
    {
      return t.m_key;
    }
  };

  // Key degree extractor.
  template <class Term>
    struct key_degree_extractor
  {
    typedef int result_type;
    int operator()(const Term &t) const
    {
      p_assert((arg_manager<Term>::assigned()));
      return t.m_key.get_degree();
    }
  };

  // Key minimum degree extractor.
  template <class Term>
    struct key_min_degree_extractor
  {
    typedef int result_type;
    int operator()(const Term &t) const
    {
      p_assert((arg_manager<Term>::assigned()));
      return t.m_key.get_min_degree();
    }
  };

  // Coefficient degree extractor.
  template <class Term>
    struct cf_degree_extractor
  {
    typedef int result_type;
    double operator()(const Term &t) const
    {
      p_assert((arg_manager<Term>::assigned()));
      return t.m_cf.get_degree();
    }
  };

  /// Index based on the degree of the coefficient.
  template <class Term>
    struct cf_degree_index
  {
    typedef boost::multi_index::indexed_by
    <
      // Unique because the composite key below ensures we cannot have duplicates.
      boost::multi_index::ordered_unique
      <
        boost::multi_index::composite_key
        <
          Term,
          cf_degree_extractor<Term>,
          key_extractor<Term>
        >,
        boost::multi_index::composite_key_compare
        <
          std::less<int>,
          std::less<typename Term::key_type>
        >
      >,
      boost::multi_index::hashed_unique<boost::multi_index::identity<Term> >
    >
    type;
  };

  /// Index based on the minimum degree of the key.
  template <class Term>
    struct key_min_degree_index
  {
    typedef boost::multi_index::indexed_by
    <
      // Unique because the composite key below ensures we cannot have duplicates.
      boost::multi_index::ordered_unique
      <
        boost::multi_index::composite_key
        <
          Term,
          key_min_degree_extractor<Term>,
          key_degree_extractor<Term>,
          key_extractor<Term>
        >,
        boost::multi_index::composite_key_compare
        <
          std::less<int>,
          std::less<int>,
          std::less<typename Term::key_type>
        >
      >,
      boost::multi_index::hashed_unique<boost::multi_index::identity<Term> >
    >
    type;
  };

  // Norm extractor. Operates on coefficient.
  template <class Term>
    struct norm_extractor
  {
    typedef double result_type;
    double operator()(const Term &t) const
    {
      p_assert((arg_manager<Term>::assigned()));
      return t.m_cf.norm(arg_manager<Term>::get());
    }
  };

  /// Norm-based index for series.
  template <class Term>
    struct norm_index
  {
    typedef boost::multi_index::indexed_by
    <
      boost::multi_index::ordered_unique
      <
        boost::multi_index::composite_key
        <
          Term,
          norm_extractor<Term>,
          key_extractor<Term>
        >,
        boost::multi_index::composite_key_compare
        <
          std::greater<double>,
          std::less<typename Term::key_type>
        >
      >,
      boost::multi_index::hashed_unique<boost::multi_index::identity<Term> >
    >
    type;
  };

  /// Hashed index on terms.
  template <class Term>
    struct keys_hash_index
  {
    typedef boost::multi_index::indexed_by
    <
      boost::multi_index::hashed_unique<boost::multi_index::identity<Term> >
    >
    type;
  };

  /// Sorted index on terms.
  template <class Term>
    struct keys_sorted_index
  {
    typedef boost::multi_index::indexed_by
    <
      boost::multi_index::ordered_unique<boost::multi_index::identity<Term> >
    >
    type;
  };
}

#endif

