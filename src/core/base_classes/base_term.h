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

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../stream_manager.h"
#include "../type_traits/eval_type.h"

#define __PIRANHA_BASE_TERM_TP_DECL class Cf, class Key, char Separator, class Derived, class Allocator
#define __PIRANHA_BASE_TERM_TP Cf,Key,Separator,Derived,Allocator

namespace piranha
{
  /// Base term class.
  /**
   * Simple composition of coefficient and key classes.
   */
  template <__PIRANHA_BASE_TERM_TP_DECL>
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
      base_term():m_cf(),m_key()
      {}
      /// Ctor from string.
      template <class ArgsTuple>
        base_term(const std::string &str, const ArgsTuple &args_tuple):m_cf(),m_key()
      {
        std::vector<std::string> vs;
        boost::split(vs,str,boost::is_any_of(std::string(separator)));
        if (vs.size() != 2)
        {
          throw bad_input();
        }
        else
        {
          m_cf = m_cf(boost::trim(vs[0]),args_tuple);
          m_key = m_key(boost::trim(vs[1]),args_tuple);
        }
      }
      /// Copy ctor.
      /**
       * Construct from base_term with different coefficient. Successful if coefficient can be converted.
       * @param[in] t base_term which will be copied.
       */
      template <class Cf2, class Derived2>
        base_term(const base_term<Cf2,Key,Separator,Derived2,Allocator> &t):m_cf(t.m_cf),m_key(t.m_key) {}
      // I/O.
      /// Print in plain format.
      template <class ArgsTuple>
        void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const
      {
        // Setup formatting.
        stream_manager::setup_print(out_stream);
        m_cf.print_plain(out_stream,args_tuple);
        out_stream << separator;
        m_key.print_plain(out_stream,args_tuple);
      }
      /// Print in latex format.
      template <class ArgsTuple>
        void print_latex(std::ostream &out_stream, const ArgsTuple &args_tuple) const
      {
        // Setup formatting
        stream_manager::setup_print(out_stream);
// TODO: redo this to work in a more general way.
//         m_cf.print_latex(out_stream,args_tuple);
//         out_stream << "&";
//         out_stream << "$" << m_key.phase(tv) << "$" << "&" << "$" << m_key.freq(tv) << "$" << "&";
//         m_key.print_latex(out_stream,args_tuple);
      }
      /// Print to stream.
      /**
       * Print format is set in piranha::stream_manager.
       */
      template <class ArgsTuple>
        void print(std::ostream &out_stream, const ArgsTuple &args_tuple) const
      {
        switch (stream_manager::format())
        {
          case stream_manager::plain:
            print_plain(out_stream,args_tuple);
            break;
          case stream_manager::latex:
            print_latex(out_stream,args_tuple);
        }
      }
      /// Print to screen.
      template <class ArgsTuple>
      void dump(const ArgsTuple &args_tuple) const
      {
        print(std::cout,args_tuple);
      }
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
      /// Assignment
      void assign(const base_term &t2)
      {
        if (this != &t2)
        {
          m_cf = t2.m_cf;
          m_key = t2.m_key;
        }
        return *this;
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
      /// Coefficient.
      /**
       * Marked as mutable for speedy operations under certain time-critical operations in hashed containers.
       * PLEASE NOTE: do _not_ abuse mutability.
       */
      mutable cf_type   m_cf;
      /// Key.
      key_type          m_key;
      /// Rebound allocator for term type.
      static allocator_type allocator;
      /// Separator between coefficient and key in I/O.
      static const char separator = Separator;
  };

  // Static members initializations.
  template <__PIRANHA_BASE_TERM_TP_DECL>
    typename base_term<__PIRANHA_BASE_TERM_TP>::allocator_type
    base_term<__PIRANHA_BASE_TERM_TP>::allocator;

  template <__PIRANHA_BASE_TERM_TP_DECL>
    const char base_term<__PIRANHA_BASE_TERM_TP>::separator;

  /// Overload of hash_value function for piranha::base_term.
  /**
   * The key's hash_value() method is used to calculate the term's hash value.
   */
  template <__PIRANHA_BASE_TERM_TP_DECL>
    inline size_t hash_value(const base_term<__PIRANHA_BASE_TERM_TP> &t)
  {
    return t.m_key.hash_value();
  }
}

#undef __PIRANHA_BASE_TERM_TP_DECL
#undef __PIRANHA_BASE_TERM_TP

#endif
