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

#ifndef PIRANHA_BASE_PSERIES_MANIP_H
#define PIRANHA_BASE_PSERIES_MANIP_H

#include "../../arg_manager.h"
#include "../../common_typedefs.h"                // For layout.
#include "../../config.h"                         // For (un)likely().
#include "../../p_exceptions.h"
#include "../../utils.h"                          // For class_converter.
#include "base_pseries_ta_macros.h"

namespace piranha
{
  /// Add phase to a term and insert the resulting terms in an external series.
  /**
   * A source term is added a phase, thus generating two terms without phases by trigonometric
   * addition formulas. The two resulting terms are written to a temporary term and are inserted into
   * an external series.
   * @param[in] phase numerical phase.
   * @param[in] term_type source term.
   * @param[in] tmp_term term used for temporary memorization of generated terms.
   * @param[out] retps series to which the terms will be added.
   * @see base_pseries::add_phase_to_term(const double &, iterator,term_type &, base_pseries &).
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::add_phase_to_term(const double &phase, const term_type &src,
    term_type &tmp_term, base_pseries &retps) const
  {
    tmp_term=src;
    // Store coefficient for later.
    cf_type tmp_c=src.cf();
    // Insert first term.
    tmp_term.cf().mult_by(std::cos(phase));
    it_s_index it_hint = retps.g_s_index().end();
    it_hint = retps.insert(tmp_term,it_hint);
    // Second term: change flavour and sign.
    switch (src.trig().flavour())
    {
      case true:
        tmp_term.trig().flavour()=false;
        tmp_term.cf()=tmp_c;
        tmp_term.cf().mult_by(-std::sin(phase));
        break;
      case false:
        tmp_term.trig().flavour()=true;
        tmp_term.cf()=tmp_c;
        tmp_term.cf().mult_by(std::sin(phase));
    }
    retps.insert(tmp_term,it_hint);
  }

  /// Add phase to a term and insert the resulting terms in an external series.
  /**
   * Analogous to base_pseries::add_phase_to_term(const double &, const term_type &,
   * term_type &, base_pseries &), the only difference is that it accepts an iterator as source term.
   * @see base_pseries::add_phase_to_term(const double &, const term_type &, term_type &, base_pseries &).
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::add_phase_to_term(const double &phase, iterator it,
    term_type &tmp_term, base_pseries &retps) const
  {
    add_phase_to_term(phase,*it,tmp_term,retps);
  }

  /// Insert numerical phases into terms.
  /**
   * Phases \f$ \phi_i \f$ are taken from a piranha::phase_list, which also specifies what type of
   * insertion is to be performed.
   * Phases are inserted term by term until the end of the series. If the phase list is shorter, addition
   * of phases stops when there are no more phases, if the phase list is longer phases in excess are
   * ignored.
   * @param[in] pl piranha::phase_list to be inserted.
   * @see piranha::phase_list.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::insert_phases(const phase_list &pl)
  {
    Derived tmp_ps;
    tmp_ps.merge_args(*this);
    phase_list::const_iterator it2=pl.begin();
    term_type tmp_term;
    for (iterator it1=begin();it1!=end();++it1)
    {
      if (it2!=pl.end())
      {
        switch (pl.operation())
        {
          case phase_list::add:
            add_phase_to_term(*it2,it1,tmp_term,tmp_ps);
            break;
          default:
            add_phase_to_term(*it2-it1->trig().phase(arguments().template get<1>()),it1,tmp_term,tmp_ps);
        }
        ++it2;
      }
      else
        // If there are no more phases just insert vanilla terms. Don't need the check, since coefficient
        // has not changed.

      {
        tmp_ps.insert(*it1);
      }
    }
    swap(tmp_ps);
  }

  /// Append coefficient arguments.
  /**
   * Calls base_pseries::append_args with Type = psymbol::cf.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::append_cf_arg(const psym_p p)
  {
    append_arg<psymbol::cf>(p);
  }

  /// Append trigonometric arguments.
  /**
   * Calls base_pseries::append_args with Type = psymbol::trig.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::append_trig_arg(const psym_p p)
  {
    append_arg<psymbol::trig>(p);
  }

  /// Swap the content of with another series.
  /**
   * All data members get swapped. This is an O(1) time operation.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::swap(Derived &ps2)
  {
    // Swap sets' contents
    s_series_set()->swap(*ps2.s_series_set());
    // Swap other members
    arguments().template get<0>().swap(ps2.arguments().template get<0>());
    arguments().template get<1>().swap(ps2.arguments().template get<1>());
    lin_args().swap(ps2.lin_args());
    static_cast<Derived *>(this)->swap_hook(ps2);
  }

  // --------------

  /// Merge arguments.
  /** Merge arguments with those of ps2. After the operation the number ofr arguments will
   * be equal or greater than ps2's. We need the template U because we want to be able to merge
   * symbols also with complex counterparts.
   * @param[in] ps2 piranha::base_pseries arguments are to be merged with.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::merge_args(const Derived2 &ps2)
  {
    if (static_cast<void *>(this) == static_cast<void const *>(&ps2))
    {
      std::cout << "Trying to merge with self, returning." << std::endl;
      return;
    }
    if (unlikely(!is_args_compatible(ps2)))
    {
      merge_incompatible_args(ps2);
    }
  }

  /// Merge argument sets which are known to be incompatible.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    void base_pseries<__PIRANHA_BASE_PS_TP>::merge_incompatible_args(const Derived2 &ps2)
  {
    Derived retval;
    // Assign old arguments and lin_args.
    retval.arguments()=arguments();
    retval.lin_args()=lin_args();
    // Get new layouts.
    // TODO: layout tuples here for more genericity.
    layout_type l_cf = utils::get_layout(cf_args(),ps2.cf_args());
    layout_type l_trig = utils::get_layout(trig_args(),ps2.trig_args());
    if (l_cf.size()> cf_type::max_size or l_trig.size()> trig_type::max_size)
    {
      throw exceptions::add_arguments("Cannot apply layout, max arguments size reached.");
    }
    // Apply layouts to arguments.
    utils::apply_layout(l_cf,retval.cf_args(),ps2.cf_args());
    utils::apply_layout(l_trig,retval.trig_args(),ps2.trig_args());
    // Apply layouts to all terms.
    const iterator it_f = end();
    term_type tmp;
    it_s_index it_hint = retval.g_s_index().end();
    for (iterator it=begin();it != it_f;++it)
    {
      tmp = *it;
      tmp.cf().apply_layout(l_cf);
      tmp.trig().apply_layout(l_trig);
      // Use safe insert because when reorganizing the args layout maybe a negative multiplier
      // may have been placed in the first position.
      it_hint = retval.insert(tmp,it_hint);
    }
    // Take care of lin_args.
    utils::apply_layout(l_trig,retval.lin_args());
    // Swap content.
    swap(retval);
  }

  // TODO: move this into fourier toolbox?
  /// Cumulative crop.
  /**
   * Crop the series from the bottom so that the sum of the norms of the cropped terms is not
   * greater than delta.
   */
//   template <__PIRANHA_BASE_PS_TP_DECL>
//     inline void base_pseries<__PIRANHA_BASE_PS_TP>::cumulative_crop(const double &delta)
//   {
//     // Won't crop a nil series
//     if (length() == 0)
//     {
//       return;
//     }
//     double part_norm=0.;
//     it_s_index it=boost::prior(g_s_index().end());
//     while (1)
//     {
//       part_norm+=it->cf().norm(arguments().template get<0>());
//       if (part_norm >= delta)
//       {
//         break;
//       }
//       if (it == g_s_index().begin())
//       {
//         term_erase(it);
//         break;
//       }
//       else
//       {
//         term_erase(it);
//         --it;
//       }
//     }
//   }

  /// Crop all terms whose norm is <= delta.
//   template <__PIRANHA_BASE_PS_TP_DECL>
//     inline void base_pseries<__PIRANHA_BASE_PS_TP>::crop(const double &delta)
//   {
//     // Won't crop a nil series
//     // TODO: replace with empty()?
//     if (length() == 0)
//     {
//       return;
//     }
//     it_s_index it=boost::prior(g_s_index().end());
//     while (1)
//     {
//       if (it->cf().norm(arguments().template get<0>()) >= delta)
//       {
//         break;
//       }
//       if (it == g_s_index().begin())
//       {
//         term_erase(it);
//         break;
//       }
//       else
//       {
//         term_erase(it);
//         --it;
//       }
//     }
//   }

  /// Crop all terms from end to it_f (included).
  // Beware that no checks are made on the validity of the iterator!
//   template <__PIRANHA_BASE_PS_TP_DECL>
//     inline void base_pseries<__PIRANHA_BASE_PS_TP>::crop(const it_s_index &it_f)
//   {
//     // Won't crop a nil series or if the crop limit is the end of the series
//     if (length()==0 || it_f==g_s_index().end())
//     {
//       return;
//     }
//     it_s_index it=boost::prior(g_s_index().end());
//     while (1)
//     {
//       term_erase(it);
//       --it;
//       if (it==it_f)
//       {
//         // Last term to be erased
//         term_erase(it);
//         break;
//       }
//     }
//   }

  /// Crop the series to a certain precision in the spectral domain.
  /**
   * The crop limit is determined by two parameters, i.e. the achieved precision in the time domain
   * and the relative desired precision in the spectral domain (relative because it is weighted against
   * the series' norm). More details in piranha::base_pseries::sdp_cutoff.
   * @param[in] achieved_tdp achieved precision in the time domain.
   * @param[in] desired relative precision in the spectral domain.
   * @see base_pseries::sdp_cutoff for a description of the cutoff procedure.
   */
//   template <__PIRANHA_BASE_PS_TP_DECL>
//     inline void base_pseries<__PIRANHA_BASE_PS_TP>::spectral_cutoff(const double &achieved_tdp,
//     const double &desired_sdp)
//   {
//     crop(sdp_cutoff(achieved_tdp,desired_sdp));
//   }
}
#endif
