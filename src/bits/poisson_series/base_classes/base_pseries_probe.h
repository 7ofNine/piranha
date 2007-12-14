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

#ifndef PIRANHA_BASE_PSERIES_PROBE_H
#define PIRANHA_BASE_PSERIES_PROBE_H

#include <algorithm>
#include <iterator>
#include <set>

#include "../../trig_evaluator.h"
#include "base_pseries_ta_macros.h"

namespace piranha
{
/// Evaluate numerically a series at a specified time. Brute force version.
/**
 * Series is evaluated using the information about the arguments. Coefficient and trigonometric
 * part methods are called internally. This version is slow and does not cache any result of evaluation.
 * Useful for debugging purposes.
 * @param[in] value, double for the time of evaluation.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    typename base_pseries<__PIRANHA_BASE_PS_TP>::eval_type base_pseries<__PIRANHA_BASE_PS_TP>::t_eval_brute(const
    double &value) const
  {
    eval_type retval(0.);
    const size_t w=trig_width();
// Terms - start from the smallest, so that we keep good precision in the summation
    const r_it_s_index it_f=g_s_index().rend();
    for (r_it_s_index it=g_s_index().rbegin();it!=it_f;++it)
    {
      retval+=it->t_eval(value,cf_s_vec_,trig_s_vec_);
    }
// Linear arguments
    for (size_t j=0;j<w;++j)
    {
      retval+=lin_args()[j]*trig_s_vec_[j]->t_eval(value);
    }
    return retval;
  };

/// Evaluate numerically a series at a specified time. Smarter version.
/**
 * Similar to base_pseries::t_eval_brute, with the different that this kind of evaluation caches evaluation of complex
 * exponential of arguments. Hence it is faster.
 * @param[in] value, double for the time of evaluation.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    typename base_pseries<__PIRANHA_BASE_PS_TP>::eval_type base_pseries<__PIRANHA_BASE_PS_TP>::t_eval(const
    double &value) const
  {
    trig_evaluator<base_pseries> te(this,value);
    eval_type retval(0.);
// Terms - start from the smallest, so that we keep good precision in the summation
    const r_it_s_index it_f=g_s_index().rend();
    for (r_it_s_index it=g_s_index().rbegin();it!=it_f;++it)
    {
      retval+=it->t_eval(te);
    }
// Linear arguments
    const size_t w=trig_width();
    for (size_t j=0;j<w;++j)
    {
      retval+=lin_args()[j]*trig_s_vec_[j]->t_eval(value);
    }
    return retval;
  };

/// Compatibility check for arguments.
/**
 * Test whether series' arguments are compatible with those from ps2. Compatibility
 * means that the first n arguments are equal, where n is the number of arguments of the series with
 * fewer arguments.
 * @param[in] ps2 piranha::base_pseries compatibility is tested against.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline bool base_pseries<__PIRANHA_BASE_PS_TP>::is_compatible(const Derived2 &ps2) const
  {
    size_t minwidth=math::min(cf_width(),ps2.cf_width()), j;
    for (j=0;j<minwidth;++j)
    {
      if (!(cf_s_vec_[j]==ps2.cf_s_vec()[j]))
      {
        return false;
      }
    }
    minwidth=math::min(trig_width(),ps2.trig_width());
    for (j=0;j<minwidth;++j)
    {
      if (!(trig_s_vec_[j]==ps2.trig_s_vec()[j]))
      {
        return false;
      }
    }
    return true;
  }

/// Check whether arguments are different.
/**
 * Returns true if the intersection of the two series' sets of arguments (both coefficient and trigonometric)
 * is void.
 */
// NOTICE: not inlined, this should not be called often and hence it would just end up increasing binary size.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    bool base_pseries<__PIRANHA_BASE_PS_TP>::args_different(const Derived2 &ps2) const
  {
// Even if there may be duplicate arguments in cf/trig_s_vec_, we don't want to use multiset:
// we are only interested in the presence or not of that argument. If there are duplicate arguments
// the insertion functions below will simply fail silently.
    typedef std::set
      <psym_p,psym_p_cmp> psym_set;
    psym_set set1, set2, set_res;
// Populate the two sets of pointer to psymbols.
    size_t j, w;
    w=cf_width();
    for (j=0;j<w;++j)
    {
      set1.insert(cf_s_vec_[j]);
    }
    w=trig_width();
    for (j=0;j<w;++j)
    {
      set1.insert(trig_s_vec_[j]);
    }
    w=ps2.cf_width();
    for (j=0;j<w;++j)
    {
      set2.insert(ps2.cf_s_vec()[j]);
    }
    w=ps2.trig_width();
    for (j=0;j<w;++j)
    {
      set2.insert(ps2.trig_s_vec()[j]);
    }
// Compute intersection of the two sets of psymbol pointers.
    std::set_intersection(set1.begin(),set1.end(),
      set2.begin(),set2.end(),
      std::inserter(set_res,set_res.begin()),
      psym_p_cmp());
    return set_res.empty();
  }

/// Calculate and return norm.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline double base_pseries<__PIRANHA_BASE_PS_TP>::g_norm() const
  {
    double retval=0.;
    const it_h_index it_f=g_h_index().end();
    for (it_h_index it=g_h_index().begin();it!=it_f;++it)
    {
      retval+=it->g_cf()->norm(cf_s_vec_);
    }
    return retval;
  }

// Return an iterator pointing to the last term before the worst discontinuity in the series
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index base_pseries<__PIRANHA_BASE_PS_TP>::discontinuity() const
  {
// We need at least 3 elements to detect a discontinuity
    if (length()<3)
    {
      return g_s_index().end();
    }
    const it_s_index it_f=boost::prior(g_s_index().end());
    it_s_index it, it_candidate=g_s_index().end();
    double rel_delta, candidate=0.;
    size_t index=0;
    for (it=g_s_index().begin();it!=it_f;++it)
    {
      rel_delta=(it->g_cf()->norm(cf_s_vec_)-boost::next(it)->g_cf()->norm(cf_s_vec_))/it->g_cf()->norm(cf_s_vec_);
      if (rel_delta>candidate)
      {
        std::cout << "Found discontinuity candidate at index position: " <<
          index << std::endl;
        candidate=rel_delta;
        it_candidate=it;
      }
      ++index;
    }
    return it_candidate;
  }

/// Establish the cutoff level for the desired precision in the spectral domain.
/** A time comparison (tc) tells us what is the maximum error against the "exact result" over a
 * certain time span. We assume that this error is split equally over all terms of a series.
 * This is called <EM>STE</EM>, or single-term error.
 *
 * A tc can be used to truncate a series (at the very end of manipulations) to cut off
 * terms based on the kind of precision that we want to achieve on their amplitudes (spectral
 * domain precision - <EM>SDP</EM>). As we move towards the end of the series the sdp will become smaller
 * and smaller. We want to be able to ditch those term whose amplitude is a certain multiplier
 * of the ste, i.e. the ste must be a certain fraction of the term's amplitude.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index base_pseries<__PIRANHA_BASE_PS_TP>::sdp_cutoff(
    const double &achieved_tdp_, const double &desired_sdp_) const
  {
    if (length()==0)
    {
      return g_s_index().end();
    }
// Take absolute values, just to prevent disasters in case of invalid input
    const double desired_sdp=std::abs(desired_sdp_), achieved_tdp=std::abs(achieved_tdp_);
// We assume each term gets an equal share of the blame :)
    const double ste=achieved_tdp/length();
    std::cout << "STE is " << ste << std::endl;
    it_s_index it;
    for (it=g_s_index().begin();it!=g_s_index().end();++it)
    {
      if (it->g_cf()->norm(cf_s_vec_)*desired_sdp < ste)
      {
        break;
      }
    }
    return it;
  }

/// Find index of argument by its name.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline size_t base_pseries<__PIRANHA_BASE_PS_TP>::trig_index(const std::string &name) const
  {
    size_t i;
    for (i=0;i<trig_width();++i)
    {
      if (trig_s_vec_[i]->name()==name)
      {
        break;
      }
    }
    if (i==trig_width()+1)
    {
      std::cout << "Error: no symbol with that name" << std::endl;
    }
    return i;
  }

/// Find the mean value of a series' evaluation over a timespan.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::eval_type base_pseries<__PIRANHA_BASE_PS_TP>::mean(
    const double &t0, const double &t1, const size_t &n) const
  {
    if (n==0)
    {
      std::cout << "Warning: won't calculate mean with zero values, returning 0." << std::endl;
      return 0.;
    }
    Derived const *derived_cast=static_cast<Derived const *>(this);
    double step=(t1-t0)/n, t=t0;
    eval_type retval(0.);
    for(size_t i=0;i<n;++i)
    {
      retval+=derived_cast->t_eval(t);
      t+=step;
    }
    return retval/(double)n;
  }

/// Memory footprint of the series.
/**
 * Internally it invokes the footprint methods of coefficients and trigonometric parts.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline size_t base_pseries<__PIRANHA_BASE_PS_TP>::footprint() const
  {
    size_t retval=sizeof(self)+trig_width()*(sizeof(psymbol *)+sizeof(int16));
    it_h_index it_f=g_h_index().end();
    for (it_h_index it=g_h_index().begin();it!=it_f;++it)
    {
      retval+=it->footprint();
    }
    return retval;
  }

/// Diagnostic check on terms.
/**
 * This functions calls Term::checkup on all terms of the series. If an error is
 * encountered it returns false, otherwise it will return true.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline bool base_pseries<__PIRANHA_BASE_PS_TP>::checkup() const
  {
    it_s_index it_f=g_s_index().end();
    for (it_s_index it=g_s_index().begin();it!=it_f;++it)
    {
      if (!it->checkup(*this))
      {
        return false;
      }
    }
    std::cout << "All Ok." << std::endl;
    return true;
  }

/// Check for single coefficient series.
/**
 * Returns true if series contains one single cosine term with null trigonometric part.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline bool base_pseries<__PIRANHA_BASE_PS_TP>::is_cf() const
  {
    if (length()==1 && g_s_index().begin()->g_trig()->g_flavour() && g_s_index().begin()->g_trig()->is_zero())
    {
      return true;
    }
    return false;
  };

/// Check whether a series is empty or not.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline bool base_pseries<__PIRANHA_BASE_PS_TP>::is_empty() const
  {
    return g_series_set()->empty();
  }

/// Determine trigonometric density.
/**
 * Trigonometric density is defined as the ratio non-zero-arguments/total-arguments.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline double base_pseries<__PIRANHA_BASE_PS_TP>::trig_density() const
  {
    const iterator it_f=end();
    double retval=0;
    size_t count=0;
    for (iterator it=begin();it!=it_f;++it)
    {
      retval+=it->g_trig()->density(*this);
      ++count;
    }
    return (retval/count);
  }
}
#endif
