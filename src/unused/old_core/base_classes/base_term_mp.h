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
  template <int LastIndex, int Index = 0> struct term_helper_brute_evaluation
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class ArgumentsTuple> static typename Term::eval_type run(const Term &term,
      const double &t, const ArgumentsTuple &a)
    {
      BOOST_STATIC_ASSERT(Term::size == (size_t)boost::tuples::length<ArgumentsTuple>::value);
      return term.elements.template get<Index>().t_eval(t, a.template get<Index>())*
        term_helper_brute_evaluation<LastIndex, Index+1>::run(term,t,a);
    }
  };

  template <int LastIndex> struct term_helper_brute_evaluation<LastIndex, LastIndex>
  {
    template <class Term, class ArgumentsTuple> static typename Term::eval_type run(const Term &term,
      const double &t, const ArgumentsTuple &a)
    {
      BOOST_STATIC_ASSERT(Term::size == (size_t)boost::tuples::length<ArgumentsTuple>::value);
      return term.elements.template get<LastIndex>().t_eval(t, a.template get<LastIndex>());
    }
  };

  template <int LastIndex> struct term_helper_checkup
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class Series> static bool run(const Term &term, const Series &s)
    {
      return (term.elements.template get<LastIndex>().checkup(s)
        and term_helper_checkup<LastIndex-1>::run(term,s));
    }
  };

  template <> struct term_helper_checkup<0>
  {
    template <class Term, class Series> static bool run(const Term &term, const Series &s)
    {
      return term.elements.template get<0>().checkup(s);
    }
  };

  template <int LastIndex> struct term_helper_ignorability
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class ArgsTuple> static bool run(const Term &term, const ArgsTuple &args_tuple)
    {
      return (term.elements.template get<LastIndex>().is_ignorable(args_tuple)
        or
        term_helper_ignorability<LastIndex-1>::run(term,args_tuple));
    }
  };

  template <> struct term_helper_ignorability<0>
  {
    template <class Term, class ArgsTuple> static bool run(const Term &term, const ArgsTuple &args_tuple)
    {
      return term.elements.template get<0>().is_ignorable(args_tuple);
    }
  };

  template <int LastIndex> struct term_helper_insertability
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class ArgsTuple> static bool run(const Term &term, const ArgsTuple &args_tuple)
    {
      return (term.elements.template get<LastIndex>().is_insertable(args_tuple)
        and
        term_helper_insertability<LastIndex-1>::run(term,args_tuple));
    }
  };

  template <> struct term_helper_insertability<0>
  {
    template <class Term, class ArgsTuple> static bool run(const Term &term, const ArgsTuple &args_tuple)
    {
      return term.elements.template get<0>().is_insertable(args_tuple);
    }
  };

  template <int LastIndex> struct term_helper_needs_padding
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class ArgsTuple> static bool run(const Term &term, const ArgsTuple &args_tuple)
    {
      return (term.elements.template get<LastIndex>().needs_padding(args_tuple)
        or
        term_helper_needs_padding<LastIndex-1>::run(term,args_tuple));
    }
  };

  template <> struct term_helper_needs_padding<0>
  {
    template <class Term, class ArgsTuple> static bool run(const Term &term, const ArgsTuple &args_tuple)
    {
      return term.elements.template get<0>().needs_padding(args_tuple);
    }
  };

  template <int LastIndex> struct term_helper_pad_right
  {
    BOOST_STATIC_ASSERT(LastIndex >= 0);
    template <class Term, class ArgsTuple> static void run(Term &term, const ArgsTuple &args_tuple)
    {
      term.elements.template get<LastIndex>().pad_right(args_tuple);
      term_helper_pad_right<LastIndex-1>::run(term,args_tuple);
    }
  };

  template <> struct term_helper_pad_right<0>
  {
    template <class Term, class ArgsTuple> static void run(Term &term, const ArgsTuple &args_tuple)
    {
      term.elements.template get<0>().pad_right(args_tuple);
    }
  };
}
#endif
