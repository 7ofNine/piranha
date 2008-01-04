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

#ifndef PIRANHA_BASE_TERM_H
#define PIRANHA_BASE_TERM_H

#include <boost/tuple/tuple.hpp>

#include "base_term_mp.h"

namespace piranha
{
  /// Base term class.
  /**
   * Wraps around a boost::tuple (http://www.boost.org/libs/tuple/doc/tuple_users_guide.html)
   * that will contain the elements of the term. Supports a maximum of 7 elements.
   */
  template < class T0, class T1 =boost::tuples::null_type, class T2 =boost::tuples::null_type,
    class T3 =boost::tuples::null_type, class T4 =boost::tuples::null_type,
    class T5 =boost::tuples::null_type, class T6 =boost::tuples::null_type> class base_term
  {
    public:
      /// Alias for the tuple type.
      typedef boost::tuple<T0,T1,T2,T3,T4,T5,T6> tuple_type;
      /// Alias for evaluation type.
      /**
       * Evaluation type is determined by the first element of the tuple.
       */
      typedef typename eval_type<T0>::type eval_type;
      /// Number of elements.
      static const size_t size = boost::tuples::length<tuple_type>::value;
      /// Last index of tuple.
      static const size_t last_index = size-1;
      /// Helper struct for the selection of element number N.
      template <int N> struct nth_element
      {
          // Make extra-sure we don't use an invalid index.
          BOOST_STATIC_ASSERT(N >= 0 and N < size);
          typedef typename boost::tuples::element<N,tuple_type>::type type;
      };
      /// Empty ctor.
      /**
       * Default-initializes all elements of the term.
       */
      base_term()
      {
      }
      /// Ctor from elements
      base_term(const T0 &e0, const T1 &e1 = boost::tuples::null_type(),
        const T2 &e2 = boost::tuples::null_type(), const T3 &e3 = boost::tuples::null_type(),
        const T4 &e4 = boost::tuples::null_type(), const T5 &e5 = boost::tuples::null_type(),
        const T6 &e6 = boost::tuples::null_type()) :
        elements(e0, e1, e2, e3, e4, e5, e6)
      {
      }
      /// Copy ctor.
      /**
       * Construct from base_term with different elements. Successful if elements can be converted.
       * @param[in] t base_term which will be copied.
       */
      template <class T0_2, class T1_2, class T2_2, class T3_2, class T4_2, class T5_2, class T6_2>
      base_term(const base_term<T0_2,T1_2,T2_2,T3_2,T4_2,T5_2,T6_2> &t):elements(t.elements)
      {
      }
      /// Numerical evaluation, brute force version.
      /**
       * Evaluate numerically the term given the time of evaluation and a vector of piranha::psymbol describing
       * the arguments. The evaluation is "dumb", in the sense that it happens term by term without caching and reusing
       * any previous calculation. Slow but reliable, hence useful for debugging purposes.
       * @param[in] t time of evaluation.
       * @param[in] a tuple of arguments relative to the elements of the term.
       */
      template <class ArgumentsTuple> eval_type t_eval_brute(const double &t, const ArgumentsTuple &a) const
      {
        return term_helper_brute_evaluation<last_index>::run(*this, t, a);
      }
      /// Run diagnostic test.
      /**
       * Run a check on the elements of the term. The exact nature of the check
       * depends on the Series template parameter. Returns true if none of the elements exhibits issues,
       * otherwise returns false.
       *
       * @param[in] s series used to retrieve checkup criterions from.
       */
      template <class Series> bool checkup(const Series &s) const
      {
        return term_helper_checkup<last_index>::run(*this, s);
      }
      /// See if a term can be ignored when inserting it into a series.
      /**
       * Returns true if at least one of the elements of the term is ignorable.
       *
       * @param[in] s series used to retrieve ignorability criterions from.
       */
      template <class Series> bool is_ignorable(const Series &s) const
      {
        return term_helper_ignorability<last_index>::run(*this, s);
      }
      /// Return const reference to the last element.
      const typename boost::tuples::element<boost::tuples::length<tuple_type>::value-1,tuple_type>::type & last_element() const
      {
        return elements.template get<last_index>();
      }
      /// Equality test.
      /**
       * By default equality is defined by the equality of the last elements of the term.
       */
      bool operator==(const base_term &t) const
      {
        return last_element() == t.last_element();
      }
      // Data members.
      /// Elements of the term.
      tuple_type elements;
  };

  /// Overload of hash_value function for piranha::base_term.
  /**
   * By default the last elements' hash_value() method is used to calculate the term's hash value.
   */
  template <class T0, class T1, class T2, class T3, class T4, class T5, class T6> inline size_t hash_value(
    const base_term<T0,T1,T2,T3,T4,T5,T6> &t)
  {
    return t.last_element().hash_value();
  }
}

#endif
