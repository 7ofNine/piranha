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

#ifndef PIRANHA_SIMPLE_TERM_H
#define PIRANHA_SIMPLE_TERM_H

#include "../../base_classes/base_term.h"

namespace piranha
{
  /// Simple Poisson series term class.
  template <class Cf, class Trig, class Allocator>
    class simple_term: public base_term<Cf,Trig,Allocator,simple_term<Cf,Trig,Allocator> >
  {
    public:
      /// Alias for the ancestor.
      typedef base_term<Cf,Trig,Allocator,simple_term<Cf,Trig,Allocator> > ancestor;
      /// Alias for coefficient type.
      typedef Cf cf_type;
      /// Alias for trigonometric type.
      typedef Trig trig_type;
      /// Alias for evaluation type.
      typedef typename ancestor::eval_type eval_type;
      /// Default constructor.
      explicit simple_term():
      ancestor::base_term() {}
      /// Constructor from generic coefficient and fixed trigonometric part.
      /**
       * Constructs from generic coefficient type.
       */
      template <class Cf2>
        explicit simple_term(const Cf2 &c, const trig_type &t):
        ancestor(cf_type(c),t)
        {}
      /// Generic copy constructor.
      /**
       * Constructs from piranha::simple_term with optionally different coefficient type.
       */
      template <class Cf2>
        explicit simple_term(const simple_term<Cf2,Trig,Allocator> &term):
        ancestor(term)
        {}
      // Getters
      /// Get mutable coefficient reference.
      cf_type &cf()
      {
        return ancestor::m_cf;
      }
      /// Get const reference to coefficient.
      const cf_type &cf() const
      {
        return ancestor::m_cf;
      }
      /// Get mutable reference to trigonometric part.
      trig_type &trig()
      {
        return ancestor::m_key;
      }
      /// Get const reference to trigonometric part.
      const trig_type &trig() const
      {
        return ancestor::m_key;
      }
      size_t footprint() const
      {
        return (sizeof(simple_term));
      }
      /// Assignment operator.
      simple_term &operator=(const simple_term &t2)
      {
        if (this != &t2)
        {
          ancestor::m_cf = t2.m_cf;
          ancestor::m_key = t2.m_key;
        }
        return *this;
      }
      /// Invert the sign of trigonometric multipliers.
      void invert_trig_args()
      {
        trig().invert_sign();
        if (!(trig().flavour()))
        {
          cf().invert_sign();
        }
      }
      // I/O.
      /// Print in plain format.
      void print_plain(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        // Setup formatting.
        stream_manager::setup_print(out_stream);
        cf().print_plain(out_stream, cv);
        out_stream << stream_manager::data_separator();
        trig().print_plain(out_stream, tv);
      }
      /// Print in latex format.
      void print_latex(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        // Setup formatting
        stream_manager::setup_print(out_stream);
        cf().print_latex(out_stream, cv);
        out_stream << "&";
        out_stream << "$" << trig().phase(tv) << "$" << "&" << "$" << trig().freq(tv) << "$" << "&";
        trig().print_latex(out_stream, tv);
      }
      /// Print to stream.
      /**
       * Print format is set in piranha::stream_manager.
       */
      void print(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        switch (stream_manager::format())
        {
          case stream_manager::plain:
            print_plain(out_stream, cv, tv);
            break;
          case stream_manager::latex:
            print_latex(out_stream, cv, tv);
        }
      }
      /// Print to screen.
      void put(const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        print(std::cout, cv, tv);
      }
      /// Smarter numerical evaluation
      /**
       * Similar to brute force evaluation, with the difference that sine and cosine of trigonometric arguments are cached
       * and re-used over the evaluation of the series. Typically faster by a factor of 2-3, depending on the series' characteristics.
       * @param[in] te piranha::trig_evaluator object that caches complex exponentials of trigonometric arguments.
       */
      template <class TrigEvaluator> eval_type t_eval(TrigEvaluator &te) const
      {
        eval_type retval=cf().t_eval(te.value(), te.ps()->arguments().template get<0>());
        retval*=trig().t_eval(te);
        return retval;
      }
      simple_term *canonicalise(simple_term *out) const
      {
// TODO: use switch?
        if (trig().sign() < 0)
        {
          if (out == 0)
          {
            out=ancestor::term_allocator.allocate(1);
            ancestor::term_allocator.construct(out,*this);
          }
          out->invert_trig_args();
        }
        return out;
      }
  };
}
#endif
