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

#ifndef PIRANHA_BASE_PSERIES_MANIP_H
#define PIRANHA_BASE_PSERIES_MANIP_H

#include "../../arg_manager.h"
#include "../../common_typedefs.h" // For layout.
#include "../../config.h" // For (un)likely().
#include "../../stats.h"
#include "../../p_exceptions.h"
#include "base_pseries_ta_macros.h"

namespace piranha
{
/// Find term.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_h_index
    base_pseries<__PIRANHA_BASE_PS_TP>::find_term(const term_type &t) const
  {
    return g_h_index().find(t.trig());
  }

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
    retps.insert_check_positive(tmp_term);
// Second term: change flavour and sign.
    switch (src.trig().g_flavour())
    {
      case true:
        tmp_term.trig().s_flavour()=false;
        tmp_term.cf()=tmp_c;
        tmp_term.cf().mult_by(-std::sin(phase));
        break;
      case false:
        tmp_term.trig().s_flavour()=true;
        tmp_term.cf()=tmp_c;
        tmp_term.cf().mult_by(std::sin(phase));
    }
    retps.insert_check_positive(tmp_term);
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

/// Append arguments.
/**
 * The arguments will be appended at the end of the corresponding argument vector.
 * @param[in] N arguments type.
 * @param[in] v vector_psym_p to be appended.
 * @see base_pseries::m_arguments tuple of arguments.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <psymbol::type Type>
    void base_pseries<__PIRANHA_BASE_PS_TP>::append_args(const vector_psym_p &v)
  {
    try
    {
      Derived retval;
      retval.lin_args()=lin_args();
      retval.arguments()=arguments();
// Append psymbols from v.
      switch (Type)
      {
// TODO: use tuple in terms here too: term_type::limits.get<Type>().max_size or similar.
        case (psymbol::cf):
          if (cf_type::max_size < v.size() + retval.arguments().template get<psymbol::cf>().size())
          {
            throw exceptions::add_arguments("Cannot append further cf_args, max_size reached.");
          }
          break;
        case (psymbol::trig):
          if (trig_type::max_size < v.size() + retval.arguments().template get<psymbol::trig>().size())
          {
            throw exceptions::add_arguments("Cannot append further trig_args, max_size reached.");
          }
// If we are dealing with trig, we must take care of lin_args too.
          retval.lin_args().insert(retval.lin_args().end(),v.size(),0);
          break;
        default:
          p_assert(false);
          ;
      }
      retval.arguments().template get<Type>().insert(retval.arguments().template get<Type>().end(),v.begin(),v.end());
// NOTICE: this can be probably split in 2 here if we want to use it in generic_series routines.
      const it_h_index it_f=g_h_index().end();
      const size_t new_size=retval.arguments().template get<Type>().size();
      for (it_h_index it=g_h_index().begin();it!=it_f;++it)
      {
// NOTICE: find a way to avoid resizes here?
        term_type tmp_term=(*it);
        switch (Type)
        {
// TODO: This would be the ideal place to use tuples inside terms...
          case (psymbol::cf):
            tmp_term.cf().pad_right(new_size);
            break;
          case (psymbol::trig):
            tmp_term.trig().pad_right(new_size);
            break;
          default:
            p_assert(false);
            ;
        }
// NOTICE: use hinted insertion here?
        retval.insert_check_positive(tmp_term);
      }
      swap(retval);
    }
    catch (exceptions::add_arguments &e) {}
  }

/// Append coefficient arguments.
/**
 * Calls base_pseries::append_args with Type = psymbol::cf.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::append_cf_args(const vector_psym_p &v)
  {
    append_args<psymbol::cf>(v);
  }

/// Append trigonometric arguments.
/**
 * Calls base_pseries::append_args with Type = psymbol::trig.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::append_trig_args(const vector_psym_p &v)
  {
    append_args<psymbol::trig>(v);
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
    lin_args_.swap(ps2.lin_args_);
    static_cast<Derived *>(this)->swap_hook(ps2);
  }


// Insert a new term into the series
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <bool Sign>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index
    base_pseries<Cf, Trig, Term, I, Derived, Allocator>::term_insert_new(const term_type &term,
    const it_s_index *it_hint)
  {
    arg_manager::arg_assigner aa(&arguments().template get<0>(),&arguments().template get<1>());
    it_s_index it_new;
    if (it_hint==0)
    {
      std::pair<iterator,bool> result=s_s_index().insert(term);
      p_assert(result.second);
      it_new=result.first;
    }
    else
    {
// FIXME: use asserts here? The problem here is that we are using hinted
// insertion, the return value is different from above (but above an assert
// is needed too).
      it_new=s_s_index().insert(*it_hint,term);
      p_assert(it_new!=end());
    }
    if (!Sign)
    {
// This is an O(1) operation, since the order in the set is not changed
// There is a re-hash involved, it still should be cheaper than
// creating a new term though.
      cf_type new_c=it_new->cf();
      new_c.invert_sign();
      action_assert(s_s_index().modify(it_new,modifier_update_cf(new_c)));
    }
    static_cast<Derived *>(this)->new_term_post_insertion_hook(term);
    return it_new;
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::term_erase(const it_h_index &it)
  {
    arg_manager::arg_assigner aa(&arguments().template get<0>(),&arguments().template get<1>());
    static_cast<Derived *>(this)->term_pre_erase_hook(*it);
    s_h_index().erase(it);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::term_erase(const it_s_index &it)
  {
    arg_manager::arg_assigner aa(&arguments().template get<0>(),&arguments().template get<1>());
    static_cast<Derived *>(this)->term_pre_erase_hook(*it);
    s_s_index().erase(it);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::term_update(const it_h_index &it, cf_type &new_c)
  {
    arg_manager::arg_assigner aa(&arguments().template get<0>(),&arguments().template get<1>());
    static_cast<Derived *>(this)->term_pre_update_hook(*it,new_c);
// Update the existing term
    action_assert(s_h_index().modify(it,modifier_update_cf(new_c)));
  }

// **************** //
// INSERT FUNCTIONS //
// **************** //
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <bool Sign>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index base_pseries<__PIRANHA_BASE_PS_TP>::ll_insert(
    const term_type &term, const it_s_index *it_hint)
  {
// We must check for ignorability here, before assertions, since at this point cf_width could be zero
// for polynomials, and thus assertion would wrongly fail.
    if (term.is_ignorable(*this))
    {
      return g_s_index().end();
    }
    p_assert(term.cf().is_insertable(cf_width()) and !term.cf().needs_padding(cf_width()));
    p_assert(term.trig().is_insertable(trig_width()) and !term.trig().needs_padding(trig_width()));
    p_assert(term.trig().sign()>0);
    it_s_index ret_it;
    it_h_index it(find_term(term));
    if (it == g_h_index().end())
    {
// The term is NOT a duplicate, insert in the set. Record where we inserted,
// so it can be used in additions and multiplications.
      ret_it=term_insert_new<Sign>(term,it_hint);
      stats::insert();
    }
    else
    {
// The term is in the set, hence an existing term will be modified.
// Add or subtract according to request.
      cf_type new_c;
      if (Sign)
      {
        new_c=it->cf();
        new_c.add(term.cf());
      }
      else
      {
        new_c=it->cf();
        new_c.subtract(term.cf());
      }
// Check if the resulting coefficient can be ignored (ie it is small).
      if (new_c.is_ignorable(*this))
      {
        term_erase(it);
      }
      else
      {
        term_update(it,new_c);
      }
// If we are erasing or updating there's no point in giving an hint on where
// the action took place, just return the end() iterator.
      ret_it=g_s_index().end();
      stats::pack();
    }
    return ret_it;
  }

  template <template <class, class> class Term, class Cf2, class Cf1, class Trig>
    struct transform_term
  {
    explicit transform_term(const Term<Cf2,Trig> &term):result(term.cf(),term.trig()) {}
    Term<Cf1,Trig>      result;
  };

  template <template <class, class> class Term, class Cf, class Trig>
    struct transform_term<Term,Cf,Cf,Trig>
  {
    explicit transform_term(const Term<Cf,Trig> &term):result(term) {}
    const Term<Cf,Trig> &result;
  };

/// Main insertion function.
/**
 * This function is used to insert terms into a series. It requires that the arguments
 * in the coefficient and in the trigonometric part of the term are fewer or as many as the series'
 * ones, otherwise an assertion fails and the program aborts. base_pseries::merge_args,
 * base_pseries::append_cf_args, base_pseries::append_trig_args, etc. can be used to add the needed arguments
 * to the series.
 *
 * This function performs some checks and then calls base_pseries::ll_insert.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Cf2, bool CheckTrigSign, bool Sign>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index base_pseries<__PIRANHA_BASE_PS_TP>::insert(
    const Term<Cf2, trig_type> &term_, const it_s_index *it_hint)
  {
    transform_term<Term,Cf2,cf_type,trig_type> term(term_);
    const size_t cw=cf_width(), tw=trig_width();
// It should not happen because resizing in this case should already be managed
// by external routines (merge_args, and input from file).
    p_assert(term.result.cf().is_insertable(cw));
    p_assert(term.result.trig().is_insertable(tw));
    term_type *new_term(0);
    const bool padding_needed=(term.result.cf().needs_padding(cw) or term.result.trig().needs_padding(tw));
    if (unlikely(padding_needed))
    {
      new_term=term_allocator.allocate(1);
      term_allocator.construct(new_term,term.result);
      new_term->cf().pad_right(cw);
      new_term->trig().pad_right(tw);
    }
    if (CheckTrigSign)
    {
      if (term.result.trig().sign() < 0)
      {
        if (new_term == 0)
        {
          new_term=term_allocator.allocate(1);
          term_allocator.construct(new_term,term.result);
        }
        new_term->invert_trig_args();
      }
    }
    const term_type *insert_term;
    if (new_term == 0)
    {
      insert_term=&term.result;
    }
    else
    {
      insert_term=new_term;
    }
    it_s_index ret_it=ll_insert<Sign>(*insert_term,it_hint);
    if (new_term != 0)
    {
      term_allocator.destroy(new_term);
      term_allocator.deallocate(new_term,1);
    }
    return ret_it;
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Cf2>
    inline typename base_pseries<__PIRANHA_BASE_PS_TP>::it_s_index base_pseries<__PIRANHA_BASE_PS_TP>::insert_check_positive(
    const Term<Cf2, trig_type> &term, const it_s_index *it_hint)
  {
    return insert<Cf2,true,true>(term,it_hint);
  }

// --------------

/// Merge arguments.
/** Merge arguments with those of ps2. After the operation the size of *this will
 * be equal or greater than ps2's. We need the template U because we want to be able to merge
 * symbols also with complex counterparts.
 * @param[in] ps2 piranha::base_pseries arguments are to be merged with.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::merge_args(const Derived2 &ps2)
  {
    if ((void *)this == (void *)&ps2)
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
    layout_type l_cf = utils::get_layout(cf_args(),ps2.cf_args());
    layout_type l_trig = utils::get_layout(trig_args(),ps2.trig_args());
    if (l_cf.size() > cf_type::max_size or l_trig.size() > trig_type::max_size)
    {
      throw exceptions::add_arguments("Cannot apply layout, max arguments size reached.");
    }
// Apply layouts to arguments.
    utils::apply_layout(l_cf,retval.cf_args(),ps2.cf_args());
    utils::apply_layout(l_trig,retval.trig_args(),ps2.trig_args());
// Apply layouts to all terms.
    const iterator it_f = end();
    term_type tmp;
    for (iterator it=begin();it != it_f;++it)
    {
// TODO: use hinting.
      tmp = *it;
      tmp.cf().apply_layout(l_cf);
      tmp.trig().apply_layout(l_trig);
// Use safe insert because when reorganizing the args layout maybe a negative multiplier
// may have been placed in the first position.
      retval.insert_check_positive(tmp);
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
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::cumulative_crop(const double &delta)
  {
// Won't crop a nil series
    if (length() == 0)
    {
      return;
    }
    double part_norm=0.;
    it_s_index it=boost::prior(g_s_index().end());
    while (1)
    {
      part_norm+=it->cf().norm(arguments().template get<0>());
      if (part_norm >= delta)
      {
        break;
      }
      if (it == g_s_index().begin())
      {
        term_erase(it);
        break;
      }
      else
      {
        term_erase(it);
        --it;
      }
    }
  }

/// Crop all terms whose norm is <= delta.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::crop(const double &delta)
  {
// Won't crop a nil series
// TODO: replace with empty()?
    if (length() == 0)
    {
      return;
    }
    it_s_index it=boost::prior(g_s_index().end());
    while (1)
    {
      if (it->cf().norm(arguments().template get<0>()) >= delta)
      {
        break;
      }
      if (it == g_s_index().begin())
      {
        term_erase(it);
        break;
      }
      else
      {
        term_erase(it);
        --it;
      }
    }
  }

/// Crop all terms from end to it_f (included).
// Beware that no checks are made on the validity of the iterator!
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::crop(const it_s_index &it_f)
  {
// Won't crop a nil series or if the crop limit is the end of the series
    if (length()==0 || it_f==g_s_index().end())
    {
      return;
    }
    it_s_index it=boost::prior(g_s_index().end());
    while (1)
    {
      term_erase(it);
      --it;
      if (it==it_f)
      {
// Last term to be erased
        term_erase(it);
        break;
      }
    }
  }

/// Crop the series to a certain precision in the spectral domain.
/**
 * The crop limit is determined by two parameters, i.e. the achieved precision in the time domain
 * and the relative desired precision in the spectral domain (relative because it is weighted against
 * the series' norm). More details in piranha::base_pseries::sdp_cutoff.
 * @param[in] achieved_tdp achieved precision in the time domain.
 * @param[in] desired relative precision in the spectral domain.
 * @see base_pseries::sdp_cutoff for a description of the cutoff procedure.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::spectral_cutoff(const double &achieved_tdp,
    const double &desired_sdp)
  {
    crop(sdp_cutoff(achieved_tdp,desired_sdp));
  }
}
#endif
