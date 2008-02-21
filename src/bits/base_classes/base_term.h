/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redis\bute it and/or modify  *
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

#include <memory>

namespace piranha
{
  /// Base term class.
  /**
   * Simple composition of coefficient and key classes.
   */
  template <class Cf, class Key, class Allocator, class Derived>
    class base_term
  {
    public:
      /// Alias for coefficient type.
      typedef Cf cf_type;
      /// Alias for key type.
      typedef Key key_type;
      /// Alias for allocator type.
      typedef typename Allocator::template rebind<Derived>::other allocator_type;
      /// Alias for evaluation type.
      /**
       * Evaluation type is determined by the coefficient.
       */
      typedef typename eval_type<cf_type>::type eval_type;
      /// Empty ctor.
      /**
       * Default-initializes coefficient and key.
       */
      base_term() {}
      /// Ctor from elements
      base_term(const cf_type &cf, const key_type &key):m_cf(cf),m_key(key)
        {}
      /// Copy ctor.
      /**
       * Construct from base_term with different coefficient and key. Successful if elements can be converted.
       * @param[in] t base_term which will be copied.
       */
      template <class Cf2, class Key2, class Derived2>
        base_term(const base_term<Cf2,Key2,Allocator,Derived2> &t):m_cf(t.m_cf),m_key(t.m_key) {}
      /// Numerical evaluation, brute force version.
      /**
       * Evaluate numerically the term given the time of evaluation and a tuple of arguments vectors.
       * The evaluation is "dumb", in the sense that it is performed term by term without caching and reusing
       * any previous calculation. Slow but reliable, hence useful for debugging purposes.
       * @param[in] t time of evaluation.
       * @param[in] a tuple of arguments vectors relative to the elements of the term.
       */
      template <class ArgumentsTuple>
        eval_type t_eval_brute(const double &t, const ArgumentsTuple &a) const
      {
        eval_type retval(m_cf.t_eval(t,a));
        retval*=m_key.t_eval(t,a);
        return retval;
      }
      /// Run diagnostic test.
      /**
       * Run a check on the elements of the term based on a tuple of arguments vectors.
       *
       * @param[in] a tuple of arguments vectors against which checkup is performed.
       */
      template <class ArgumentsTuple>
        bool checkup(const ArgumentsTuple &a) const
      {
        return (m_cf.checkup(a) and m_key.checkup(a));
      }
      /// Check whether a term can be ignored when inserting it into a series.
      /**
       * Returns true if at least one of the elements of the term is ignorable.
       *
       * @param[in] a tuple of arguments vectors against which ignorability is tested.
       */
      template <class ArgsTuple>
        bool is_ignorable(const ArgsTuple &a) const
      {
        return (m_cf.is_ignorable(a) or m_key.is_ignorable(a));
      }
      /// Check whether a term can be inserted into a series.
      /**
       * Returns true if the number of arguments provided by the arguments tuple for each element of the term
       * is compatible for insertion of the term into a series, false otherwise.
       *
       * @param[in] a tuple of arguments vectors used for the check.
       */
      template <class ArgsTuple>
        bool is_insertable(const ArgsTuple &a) const
      {
        return (m_cf.is_insertable(a) and m_key.is_insertable(a));
      }
      /// Check whether a term needs padding for insertion into a series.
      /**
       * Returns true if term needs right-padding with zeroes before insertion into series,
       * false otherwise.
       *
       * @param[in] a tuple of arguments vectors used for the check for padding.
       */
      template <class ArgsTuple>
        bool needs_padding(const ArgsTuple &a) const
      {
        return (m_cf.needs_padding(a) or m_key.needs_padding(a));
      }
      /// Pad right all elements of the term.
      /**
       * Padding sizes are taken from arguments tuple.
       *
       * @param[in] a tuple of arguments vectors used for padding.
       */
      template <class ArgsTuple>
        void pad_right(const ArgsTuple &a)
      {
        m_cf.pad_right(a);
        m_key.pad_right(a);
      }
      /// Equality test.
      /**
       * Equality is defined by the equality of the keys.
       */
      bool operator==(const base_term &t) const
      {
        return (m_key == t.m_key);
      }
      /// Hasher functor.
      /**
       * Useful in STL-like containers.
       */
      struct hasher
      {
        size_t operator()(const base_term &t) const
        {
          return t.m_key.hash_value();
        }
      };
      // Data members.
      /// Elements of the term.
      /**
       * Marked as mutable for speedy operations under certain time-critical operations in hashed containers.
       * PLEASE NOTE: do _not_ abuse mutability.
       */
      mutable cf_type   m_cf;
      mutable key_type  m_key;
      /// Rebound allocator for term type.
      static allocator_type allocator;
  };


  template <class Cf, class Key, class Allocator, class Derived>
    typename base_term<Cf,Key,Allocator,Derived>::allocator_type
    base_term<Cf,Key,Allocator,Derived>::allocator;

  /// Overload of hash_value function for piranha::base_term.
  /**
   * The key's hash_value() method is used to calculate the term's hash value.
   */
  template <class Cf, class Key, class Allocator, class Derived>
    inline size_t hash_value(const base_term<Cf,Key,Allocator,Derived> &t)
  {
    return t.m_key.hash_value();
  }
}

#endif
