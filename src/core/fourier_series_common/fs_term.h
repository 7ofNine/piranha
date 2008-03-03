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

#ifndef PIRANHA_FS_TERM_H
#define PIRANHA_FS_TERM_H

#include <string>

#include "../base_classes/base_term.h"

namespace piranha
{
  /// Term class for Fourier series.
  template <class Cf, class Trig, char Separator, class Allocator>
    class fs_term: public base_term<Cf,Trig,Separator,fs_term<Cf,Trig,Separator,Allocator>,Allocator>
  {
      /// Alias for the ancestor.
      typedef base_term<Cf,Trig,Separator,fs_term<Cf,Trig,Separator,Allocator>,Allocator> ancestor;
      /// Alias for evaluation type.
      typedef typename ancestor::eval_type eval_type;
    public:
      /// Alias for coefficient type.
      typedef Cf cf_type;
      /// Alias for trigonometric type.
      typedef Trig trig_type;
      /// Default constructor.
      explicit fs_term():ancestor::base_term() {}
      /// Ctor from string.
      template <class ArgsTuple>
        explicit fs_term(const std::string &str, const ArgsTuple &args_tuple):
        ancestor::base_term(str,args_tuple)
      {}
      /// Constructor from generic coefficient and fixed trigonometric part.
      /**
       * Constructs from generic coefficient type.
       */
      template <class Cf2>
        explicit fs_term(const Cf2 &c, const trig_type &t):ancestor(cf_type(c),t)
      {}
      /// Generic copy constructor.
      /**
       * Constructs from piranha::fs_term with optionally different coefficient type.
       */
      template <class Cf2>
        explicit fs_term(const fs_term<Cf2,Trig,Separator,Allocator> &term):ancestor(term)
      {}
      /// Smarter numerical evaluation
      /**
       * Similar to brute force evaluation, with the difference that sine and cosine of trigonometric arguments are cached
       * and re-used over the evaluation of the series. Typically faster by a factor of 2-3, depending on the series' characteristics.
       * @param[in] te piranha::trig_evaluator object that caches complex exponentials of trigonometric arguments.
       */
      template <class TrigEvaluator>
        eval_type t_eval(TrigEvaluator &te) const
      {
        eval_type retval=ancestor::m_cf.t_eval(te.m_value,te.m_ps->m_arguments);
        retval*=ancestor::m_key.t_eval(te);
        return retval;
      }
      /// Check if the term is canonical.
      bool is_canonical() const
      {
        return (ancestor::m_key.sign() > 0);
      }
      // TODO: check if it makes sense to skip the check here and assume canonicalise will be used iff
      // is_canonical has already been tested.
      /// Canonicalise the term.
      void canonicalise()
      {
        if (!is_canonical())
        {
          invert_trig_args();
        }
      }
    private:
      // Invert the sign of trigonometric multipliers.
      void invert_trig_args()
      {
        ancestor::m_key.invert_sign();
        if (!(ancestor::m_key.flavour()))
        {
          ancestor::m_cf.invert_sign();
        }
      }
  };
}
#endif
