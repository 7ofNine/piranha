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

#ifndef PIRANHA_PS_TERM_PROBE_H
#define PIRANHA_PS_TERM_PROBE_H

#include "arg_manager.h"                          // ps_term::operator<.

namespace piranha
{
/// See if a term can be ignored when inserting it into a series.
/**
 * @param[in] v vector of piranha::psymbol objects for the coefficient.
 */
  template <class Cf,class Trig>
    inline bool ps_term<Cf,Trig>::is_ignorable(const vector_psym_p &v) const
  {
// First check: numerical coefficient
// NOTE: store and use norm here?
// NOTE: the norm should go in the coefficient manipulator...
    if (g_cf().is_zero(v))
    {
      return true;
    }
// Second check: if sine, check that there is at least 1 non-zero trig arg. Otherwise ignore.
    if (g_flavour()==false)
    {
      return g_trig().is_zero();
    }
    return false;
  }

/// Sign of first non-zero trigonometric argument.
  template <class Cf,class Trig>
    inline int ps_term<Cf,Trig>::trig_sign() const
  {
    return g_trig().sign();
  }

/// Numerical evaluation.
/**
 * Evaluate numerically the term given the time of evaluation and a vector of piranha::psymbol describing
 * the arguments.
 * @param[in] t time of evaluation.
 * @param[in] vc vector of piranha::psymbol objects for the coefficient.
 * @param[in] vt vector of piranha::psymbol objects for the trigonometric part.
 */
  template <class Cf,class Trig>
    inline typename ps_term<Cf,Trig>::cf_type::eval_type
    ps_term<Cf,Trig>::t_eval(double t,const vector_psym_p &vc, const vector_psym_p &vt) const
  {
    typename cf_type::eval_type retval=g_cf().t_eval(t,vc);
    switch (g_flavour())
    {
      case true:
        retval*=std::cos(g_trig().t_eval(t,vt));
        break;
      case false:
        retval*=std::sin(g_trig().t_eval(t,vt));
    }
    return retval;
  }


// NOTICE: drop freq and phase from here, define only in trig_args and call from there?
/// Numerical phase of the term.
/**
 * Get the numerical phase of the term, given a vector of piranha::psymbol describing its arguments.
 * @param[in] v vector of piranha::psymbol objects.
 */
  template <class Cf,class Trig>
    inline double ps_term<Cf,Trig>::phase(const vector_psym_p &v) const
  {
    return g_trig().phase(v);
  }

/// Frequency of the term.
/**
 * Get the frequency of the term, given a vector of piranha::psymbol describing its arguments.
 * @param[in] v vector of piranha::psymbol objects.
 */
  template <class Cf,class Trig>
    inline double ps_term<Cf,Trig>::freq(const vector_psym_p &v) const
  {
    return g_trig().freq(v);
  }

/// Test for equality (used in the hashed index).
  template <class Cf,class Trig>
    inline bool operator==(const ps_term<Cf,Trig> &a,
    const ps_term<Cf,Trig> &b)
  {
    if (a.g_flavour()!=b.g_flavour())
    {
      return false;
    }
    return a.g_trig()==b.g_trig();
  }

/// Memory footprint.
  template <class Cf,class Trig>
    inline size_t ps_term<Cf,Trig>::footprint() const
  {
    return (sizeof(self)+g_trig().data_footprint());
  }

/// Diagnose problems.
/**
 * Run a check on the coefficient and on the trigonometric part. The exact nature of the check
 * depends on the template parameters implementations.
 */
  template <class Cf,class Trig>
    inline bool ps_term<Cf,Trig>::checkup(const size_t &cw, const size_t &tw) const
  {
    if (!g_cf().checkup(cw))
    {
      std::cout << "Coefficient failed checkup." << std::endl;
      return false;
    }
    if (!g_trig().checkup(tw))
    {
      std::cout << "Trigonometric part failed checkup." << std::endl;
      return false;
    }
    return true;
  }
}
#endif
