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

#include "arg_manager.h" // ps_term::operator<.

namespace piranha
  {
  /// Operator less-than.
  template <class Cf,class Trig>
  inline bool ps_term<Cf,Trig>::operator<(const ps_term &t2) const
    {
      p_assert(arg_manager::assigned());
      if (norm(*arg_manager::cf_args())!=t2.norm(*arg_manager::cf_args()))
        {
          return norm(*arg_manager::cf_args())>t2.norm(*arg_manager::cf_args());
        }
      else
        {
          if (flavour()<t2.flavour())
            {
              return true;
            }
          else if (t2.flavour()<flavour())
            {
              return false;
            }
          return (trig_args()>t2.trig_args());
        }
    }


  /// See if a term can be ingored when inserting it into a series.
  /**
   * @param[in] v vector of piranha::psymbol objects for the coefficient.
   */
  template <class Cf,class Trig>
  inline bool ps_term<Cf,Trig>::is_ignorable(const vector_psym_p &v) const
    {
      // First check: numerical coefficient
      // NOTE: store and use norm here?
      // NOTE: the norm should go in the coefficient manipulator...
      if (c_.is_zero(v))
        {
          return true;
        }
      // Second check: if sine, check that there is at least 1 non-zero trig arg. Otherwise ignore.
      if (flavour_==false)
        {
          return trig_args_.is_zero();
        }
      return false;
    }


  /// Sign of first non-zero trigonometric argument.
  template <class Cf,class Trig>
  inline int ps_term<Cf,Trig>::trig_sign() const
    {
      return trig_args_.sign();
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
      typename cf_type::eval_type retval=c_.t_eval(t,vc);
      switch (flavour_)
        {
        case true:
          retval*=std::cos(trig_args_.t_eval(t,vt));
          break;
        case false:
          retval*=std::sin(trig_args_.t_eval(t,vt));
        }
      return retval;
    }


  /// Get term's norm.
  /**
   * @param[in] vc vector of piranha::psymbol objects for the coefficient.
   */
  template <class Cf,class Trig>
  inline double ps_term<Cf,Trig>::norm(const vector_psym_p &vc) const
    {
      return c_.norm(vc);
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
      return trig_args_.phase(v);
    }

  /// Frequency of the term.
  /**
   * Get the frequency of the term, given a vector of piranha::psymbol describing its arguments.
   * @param[in] v vector of piranha::psymbol objects.
   */
  template <class Cf,class Trig>
  inline double ps_term<Cf,Trig>::freq(const vector_psym_p &v) const
    {
      return trig_args_.freq(v);
    }


  /// Test for equality (used in the hashed index).
  template <class Cf,class Trig>
  inline bool operator==(const ps_term<Cf,Trig> &a,
                         const ps_term<Cf,Trig> &b)
  {
    if (a.flavour()!=b.flavour())
      {
        return false;
      }
    return a.trig_args()==b.trig_args();
  }


  /// Memory footprint.
  template <class Cf,class Trig>
  inline size_t ps_term<Cf,Trig>::footprint() const
    {
      return (sizeof(self)+trig_args_.data_footprint());
    }


  /// Diagnose problems.
  /**
   * Run a check on the coefficient and on the trigonometric part. The exact nature of the check
   * depends on the template parameters implementations.
   */
  template <class Cf,class Trig>
  inline bool ps_term<Cf,Trig>::checkup(const size_t &cw, const size_t &tw) const
    {
      if (!c_.checkup(cw))
        {
          std::cout << "Coefficient failed checkup." << std::endl;
          return false;
        }
      if (!trig_args_.checkup(tw))
        {
          std::cout << "Trigonometric part failed checkup." << std::endl;
          return false;
        }
      return true;
    }
}
#endif
