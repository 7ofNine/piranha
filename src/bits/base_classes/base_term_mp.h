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

#ifndef PIRANHA_BASE_TERM_MP_H
#define PIRANHA_BASE_TERM_MP_H

/** @file base_term_mp.h
    \brief Meta-programming for piranha::base_term.
*/

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>

namespace piranha
{
  template <int LastIndex, int Index = 0>
    struct term_helper_brute_evaluation
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class ArgumentsTuple>
      static typename Term::eval_type run(const Term &term, const double &t, const ArgumentsTuple &a)
    {
      BOOST_STATIC_ASSERT(Term::size == (size_t)boost::tuples::length<ArgumentsTuple>::value);
      return term.elements.template get<Index>().t_eval(t,a.template get<Index>())*
        term_helper_brute_evaluation<LastIndex,Index+1>::run(term,t,a);
    }
  };

  template <int LastIndex>
    struct term_helper_brute_evaluation<LastIndex,LastIndex>
  {
    template <class Term, class ArgumentsTuple>
      static typename Term::eval_type run(const Term &term, const double &t, const ArgumentsTuple &a)
    {
      BOOST_STATIC_ASSERT(Term::size == (size_t)boost::tuples::length<ArgumentsTuple>::value);
      return term.elements.template get<LastIndex>().t_eval(t,a.template get<LastIndex>());
    }
  };

  template <int LastIndex, int Index = 0>
    struct term_helper_checkup
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class Series>
      static bool run(const Term &term, const Series &s)
    {
      return (term.elements.template get<Index>().checkup(s) and term_helper_checkup<LastIndex,Index+1>::run(term,s));
    }
  };

  template <int LastIndex>
    struct term_helper_checkup<LastIndex,LastIndex>
  {
    template <class Term, class Series>
      static bool run(const Term &term, const Series &s)
    {
      return term.elements.template get<LastIndex>().checkup(s);
    }
  };

  template <int LastIndex, int Index = 0>
    struct term_helper_ignorability
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class Series>
      static bool run(const Term &term, const Series &s)
    {
      return (term.elements.template get<Index>().is_ignorable(s) or
        term_helper_ignorability<LastIndex,Index+1>::run(term,s));
    }
  };

  template <int LastIndex>
    struct term_helper_ignorability<LastIndex,LastIndex>
  {
    template <class Term, class Series>
      static bool run(const Term &term, const Series &s)
    {
      return term.elements.template get<LastIndex>().is_ignorable(s);
    }
  };
}

#endif
