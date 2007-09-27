/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#include "trig_evaluator.h"

namespace piranha
{
/// Base Poisson series term class.
/**
 * Accepts coefficient and trigonometric part as template parameters.
 */
  template <class Cf, class Trig, class Derived>
    class base_term
  {
    public:
/// Alias for self.
      typedef base_term self;
/// Alias for coefficient type.
      typedef Cf cf_type;
/// Alias for trigonometric type.
      typedef Trig trig_type;
/// Alias for evaluation type.
      typedef typename cf_type::eval_type eval_type;
// Getters
/// Get flavour.
      bool &s_flavour()
      {
        return static_cast<Derived *>(this)->s_trig()->s_flavour();
      }
      const bool &g_flavour() const
      {
        return static_cast<Derived const *>(this)->g_trig()->g_flavour();
      }
/// Diagnose problems.
/**
 * Run a check on the coefficient and on the trigonometric part. The exact nature of the check
 * depends on the template parameters.
 */
// TODO: templatize this and pass series as arguments.
      bool checkup(const size_t &cw, const size_t &tw) const
      {
        if (!static_cast<Derived const *>(this)->g_cf()->checkup(cw))
        {
          std::cout << "Coefficient failed checkup." << std::endl;
          return false;
        }
        if (!static_cast<Derived const *>(this)->g_trig()->checkup(tw))
        {
          std::cout << "Trigonometric part failed checkup." << std::endl;
          return false;
        }
        return true;
      }
/// See if a term can be ignored when inserting it into a series.
/**
 * @param[in] v vector of piranha::psymbol objects for the coefficient.
 */
      bool is_ignorable(const vector_psym_p &v) const
      {
// First check: numerical coefficient
// NOTE: store and use norm here?
// NOTE: the norm should go in the coefficient manipulator...
        if (static_cast<Derived const *>(this)->g_cf()->is_zero(v))
        {
          return true;
        }
// TODO: move this into trig_args, once we move flavour there.
// Second check: if sine, check that there is at least 1 non-zero trig arg. Otherwise ignore.
        if (g_flavour()==false)
        {
          return static_cast<Derived const *>(this)->g_trig()->is_zero();
        }
        return false;
      }
// Manipulation.
      void invert_trig_args()
      {
        static_cast<Derived *>(this)->s_trig()->invert_sign();
        if (!(static_cast<Derived const *>(this)->g_trig()->g_flavour()))
        {
// FIXME: maybe here a invert_sign function for the coefficient should be used?
          *static_cast<Derived *>(this)->s_cf()*=-1;
        }
      }
// I/O.
/// Print in plain format.
      void print_plain(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
// Setup formatting.
        stream_manager::setup_print(out_stream);
        static_cast<Derived const *>(this)->g_cf()->print_plain(out_stream,cv);
        out_stream << stream_manager::data_separator();
        static_cast<Derived const *>(this)->g_trig()->print_plain(out_stream,tv);
        switch (g_flavour())
        {
          case true:
            out_stream << "c";
            break;
          case false:
            out_stream << "s";
        }
      }
/// Print in latex format.
      void print_latex(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
// Setup formatting
        stream_manager::setup_print(out_stream);
        static_cast<Derived const *>(this)->g_cf()->print_latex(out_stream,cv);
        out_stream << "&";
        out_stream << "$" << static_cast<Derived const *>(this)->g_trig()->phase(tv) <<
          "$" << "&" << "$" << static_cast<Derived const *>(this)->g_trig()->freq(tv) << "$" << "&";
        switch (g_flavour())
        {
          case true:
            out_stream << "c&";
            break;
          case false:
            out_stream << "s&";
        }
        static_cast<Derived const *>(this)->g_trig()->print_latex(out_stream,tv);
      }
/// Print to stream.
/**
 * Print format (i.e., plain or latex) is set in piranha::stream_manager.
 */
      void print(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        switch (stream_manager::format())
        {
          case stream_manager::plain:
            print_plain(out_stream,cv,tv);
            break;
          case stream_manager::latex:
            print_latex(out_stream,cv,tv);
        }
      }
/// Print to screen.
      void put(const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        print(std::cout,cv,tv);
      }
/// Numerical evaluation, brute force version.
/**
 * Evaluate numerically the term given the time of evaluation and a vector of piranha::psymbol describing
 * the arguments. The evaluation is "dumb", in the sense that it happens term by term without caching and reusing
 * any previous calculation. Slow but reliable, hence useful for debugging purposes.
 * @param[in] t time of evaluation.
 * @param[in] vc vector of piranha::psymbol objects for the coefficient.
 * @param[in] vt vector of piranha::psymbol objects for the trigonometric part.
 */
      eval_type t_eval(const double &t, const vector_psym_p &vc, const vector_psym_p &vt) const
      {
        eval_type retval=static_cast<Derived const *>(this)->g_cf()->t_eval(t,vc);
// TODO: move this into trig_args, once we move flavour there.
        switch (g_flavour())
        {
          case true:
            retval*=std::cos(static_cast<Derived const *>(this)->g_trig()->t_eval(t,vt));
            break;
          case false:
            retval*=std::sin(static_cast<Derived const *>(this)->g_trig()->t_eval(t,vt));
        }
        return retval;
      }
/// Smarter numerical evaluation
/**
 * Similar to brute force evaluation, with the difference that sine and cosine of trigonometric arguments are cached
 * and re-used over the evaluation of the series. Typically faster by a factor of 2-3, depending on the series' characteristics.
 * @param[in] te piranha::trig_evaluator object that caches complex exponentials of trigonometric arguments.
 */
      template <class TrigEvaluator>
        eval_type t_eval(TrigEvaluator &te) const
      {
        eval_type retval=static_cast<Derived const *>(this)->g_cf()->t_eval(te.value(),te.ps()->cf_s_vec());
// TODO: move this into trig_args, once we move flavour there.
        switch (g_flavour())
        {
          case true:
            retval*=static_cast<Derived const *>(this)->g_trig()->t_eval(te).real();
            break;
          case false:
            retval*=static_cast<Derived const *>(this)->g_trig()->t_eval(te).imag();
        }
        return retval;
      }
  };
}
#endif
