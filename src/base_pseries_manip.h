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

#include "arg_manager.h"
#include "stats.h"

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
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::add_phase_to_term(const double &phase, const term_type &src,
    term_type &tmp_term, base_pseries &retps) const
  {
    tmp_term=src;
// Store coefficient for later.
    cf_type tmp_c=src.g_cf();
// Insert first term.
    tmp_term.s_cf()*=std::cos(phase);
    retps.insert(tmp_term);
// Second term: change flavour and sign.
    switch (src.g_flavour())
    {
      case true:
        tmp_term.s_flavour()=false;
        tmp_term.s_cf()=tmp_c;
        tmp_term.s_cf()*=(-std::sin(phase));
        break;
      case false:
        tmp_term.s_flavour()=true;
        tmp_term.s_cf()=tmp_c;
        tmp_term.s_cf()*=std::sin(phase);
    }
    retps.insert(tmp_term);
  }

/// Add phase to a term and insert the resulting terms in an external series.
/**
 * Analogous to base_pseries::add_phase_to_term(const double &, const term_type &,
 * term_type &, base_pseries &), the only difference is that it accepts an iterator as source term.
 * @see base_pseries::add_phase_to_term(const double &, const term_type &, term_type &, base_pseries &).
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::add_phase_to_term(const double &phase, iterator it,
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
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::insert_phases(const phase_list &pl)
  {
    base_pseries tmp_ps;
    p_assert(tmp_ps.merge_args(*this));
    phase_list::const_iterator it2=pl.begin();
    term_type tmp_term;
    for (iterator it1=begin();it1!=end();++it1)
    {
      if (it2!=pl.end())
      {
        switch (pl.operation())
        {
          case phase_list::add
            :
          add_phase_to_term(*it2,it1,tmp_term,tmp_ps);
          break;
          default:
            add_phase_to_term(*it2-it1->phase(trig_s_vec_),it1,tmp_term,tmp_ps);
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

/// Prepend coefficient arguments.
/**
 * Prepend a vector of argument pointers to the current vector of coefficient arguments.
 * @param[in] v vector_psym_p to be prepended.
 * @see base_pseries::cf_s_vec_ vector of coefficient arguments for a series.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::prepend_cf_args(const vector_psym_p &v)
  {
    base_pseries retval=base_pseries();
    retval.lin_args_=lin_args_;
    retval.cf_s_vec_=cf_s_vec_;
    retval.trig_s_vec_=trig_s_vec_;
// Prepend psymbols from v.
    retval.cf_s_vec_.insert(retval.cf_s_vec_.begin(),v.begin(),v.end());
    const it_h_index it_f=h_index().end();
    const size_t n=v.size();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
// NOTICE: find a way to avoid resizes here?
      term_type tmp_term=(*it);
      tmp_term.s_cf().prepend_args(n);
// NOTICE: use hinted insertion here?
      retval.insert(tmp_term);
    }
    swap(retval);
  }

/// Append coefficient arguments.
/**
 * The arguments will be appended at the end of the coefficient argument vector.
 * @param[in] v vector_psym_p to be appended.
 * @see base_pseries::cf_s_vec_ vector of coefficient arguments for a series.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::append_cf_args(const vector_psym_p &v)
  {
    base_pseries retval=base_pseries();
    retval.lin_args_=lin_args_;
    retval.cf_s_vec_=cf_s_vec_;
    retval.trig_s_vec_=trig_s_vec_;
// Append psymbols from v.
    retval.cf_s_vec_.insert(retval.cf_s_vec_.end(),v.begin(),v.end());
    const it_h_index it_f=h_index().end();
    const size_t n=v.size();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
// NOTICE: find a way to avoid resizes here?
      term_type tmp_term=(*it);
      tmp_term.s_cf().append_args(n);
// NOTICE: use hinted insertion here?
      retval.insert(tmp_term);
    }
    swap(retval);
  }

/// Prepend trigonometric arguments.
/**
 * Prepend a vector of argument pointers to the current vector of trigonometric arguments.
 * @param[in] v vector_psym_p to be prepended.
 * @see base_pseries::trig_s_vec_ vector of trigonometric arguments for a series.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::prepend_trig_args(const vector_psym_p &v)
  {
    base_pseries retval=base_pseries();
    retval.lin_args_=lin_args_;
    retval.cf_s_vec_=cf_s_vec_;
    retval.trig_s_vec_=trig_s_vec_;
    const size_t n=v.size();
// Prepend psymbols from v.
    retval.trig_s_vec_.insert(retval.trig_s_vec_.begin(),v.begin(),v.end());
    retval.lin_args_.insert(retval.lin_args_.begin(),n,0);
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
// NOTICE: find a way to avoid resizes here?
      term_type tmp_term=(*it);
      tmp_term.s_trig().prepend_args(n);
// NOTICE: use hinted insertion here?
      retval.insert(tmp_term);
    }
    swap(retval);
  }

/// Append trigonometric arguments.
/**
 * The argument will be appended at the end of the trigonometric argument vector.
 * @param[in] v vector_psym_p to be appended.
 * @see base_pseries::trig_s_vec_ vector of trigonometric arguments for a series.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::append_trig_args(const vector_psym_p &v)
  {
    base_pseries retval=base_pseries();
    retval.lin_args_=lin_args_;
    retval.cf_s_vec_=cf_s_vec_;
    retval.trig_s_vec_=trig_s_vec_;
    const size_t n=v.size();
// Append psymbols from v.
    retval.trig_s_vec_.insert(retval.trig_s_vec_.end(),v.begin(),v.end());
    retval.lin_args_.insert(retval.lin_args_.end(),n,0);
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      term_type tmp_term=(*it);
      tmp_term.s_trig().append_args(n);
// NOTICE: use hinted insertion here?
      retval.insert(tmp_term);
    }
    swap(retval);
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::set_flavour(bool flavour)
  {
    base_pseries retval=base_pseries();
    retval.lin_args_=lin_args_;
    retval.cf_s_vec_=cf_s_vec_;
    retval.trig_s_vec_=trig_s_vec_;
// Insert terms from this into retval.
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      term_type tmp_term=(*it);
      tmp_term.s_flavour()=flavour;
// NOTICE: use hinted insertion here?
      retval.insert(tmp_term);
    }
    swap(retval);
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::upgrade_norm(const double &new_real)
  {
    norm_+=std::abs(new_real);
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::downgrade_norm(const double &new_real)
  {
    norm_-=std::abs(new_real);
  }

/// Swap the content of with another series.
/**
 * All data members get swapped. This is an O(1) time operation.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::swap(base_pseries<Cf, Trig, Term, I, Derived> &ps2)
  {
// Swap sets' contents
    set_.swap(ps2.set_);
// Swap other members
    std::swap(norm_,ps2.norm_);
    cf_s_vec_.swap(ps2.cf_s_vec_);
    trig_s_vec_.swap(ps2.trig_s_vec_);
    lin_args_.swap(ps2.lin_args_);
  }

// Insert a term into the series

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline typename base_pseries<Cf, Trig, Term, I, Derived>::it_s_index
    base_pseries<Cf, Trig, Term, I, Derived>::term_insert_new(const term_type &term, bool sign, const
    it_s_index *it_hint)
  {
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
    it_s_index it_new;
    if (it_hint==0)
    {
      std::pair<iterator,bool> result=s_index().insert(term);
      p_assert(result.second);
      it_new=result.first;
    }
    else
    {
// FIXME: use asserts here? The problem here is that we are using hinted
// insertion, the return value is different from above (but above an assert
// is needed too
      it_new=s_index().insert(*it_hint,term);
      p_assert(it_new!=end());
    }
    if (sign==false)
    {
// This is an O(1) operation, since the order in the set is not changed
// Yes.. but there is a re-hash involved, it still should be cheaper than
// creating a new term though. Investigate.
      cf_type new_c=it_new->g_cf();
// FIXME: change sign directly here, so we do not do a copy below!
      p_assert(s_index().modify(it_new,typename Term<cf_type,trig_type>::
        modifier_update_cf(-new_c)));
    }
// No need do differentiate between + and -, will get abs()-ed anyway
    upgrade_norm(term.g_cf().norm(cf_s_vec_));
    return it_new;
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::term_erase(const it_h_index &it)
  {
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
    downgrade_norm(it->norm(cf_s_vec_));
    h_index().erase(it);
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::term_erase(const it_s_index &it)
  {
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
    downgrade_norm(it->norm(cf_s_vec_));
    s_index().erase(it);
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::term_update(const it_h_index &it, const cf_type &new_c)
  {
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
// Delete old c from norm
    downgrade_norm(it->norm(cf_s_vec_));
// Update the existing term
    p_assert(h_index().modify(it,typename Term<cf_type,trig_type>::modifier_update_cf(new_c)));
// Update the norm with the new c
    upgrade_norm(new_c.norm(cf_s_vec_));
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::term_update(const it_s_index &it, const cf_type &new_c)
  {
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
// Delete old c from norm
    downgrade_norm(it->norm());
// Update the existing term
    p_assert(s_index().modify(it,typename Term<cf_type,trig_type>::modifier_update_cf(new_c)));
// Update the norm with the new c
    upgrade_norm(new_c.norm());
  }

// **************** //
// INSERT FUNCTIONS //
// **************** //
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline typename base_pseries<Cf, Trig, Term, I, Derived>::it_s_index base_pseries<Cf, Trig, Term, I, Derived>::ll_insert(
    const term_type &term, bool sign, const it_s_index *it_hint)
  {
    p_assert(term.g_cf().compatible(cf_width()));
    p_assert(term.g_trig().compatible(trig_width()));
    p_assert(term.g_trig().sign()>0);
    it_s_index ret_it;
    it_h_index it=h_index().find(boost::make_tuple(term.g_flavour(),
      term.g_trig()));
    if (it==h_index().end())
    {
// The term is NOT a duplicate, insert in the set. Record where we inserted,
// so it can be used in additions and multiplications
      ret_it=term_insert_new(term,sign,it_hint);
      stats::insert();
    }
    else
    {
// The term is in the set, hence an existing term will be modified
// Add or subtract according to request
      cf_type new_c;
      if (sign)
      {
        new_c=it->g_cf();
        new_c+=term.g_cf();
      }
      else
      {
        new_c=it->g_cf();
        new_c-=term.g_cf();
      }
// Check if the resulting coefficient can be ignored (ie it is small).
      if (new_c.is_zero(cf_s_vec_))
      {
        term_erase(it);
      }
      else
      {
        term_update(it,new_c);
      }
// If we are erasing or updating there's no point in giving an hint on where
// the action took place, just return the end() iterator.
      ret_it=s_index().end();
      stats::pack();
    }
    return ret_it;
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2>
    inline typename base_pseries<Cf, Trig, Term, I, Derived>::it_s_index base_pseries<Cf, Trig, Term, I, Derived>::ll_insert(
    const Term<Cf2,trig_type> &term, bool sign, const it_s_index *it_hint)
  {
//BOOST_STATIC_ASSERT(sizeof(U)==0);
    return ll_insert(term_type(term),sign,it_hint);
  }

/// Main insert function.
/**
 * This function is used to insert terms into a series. It requires that the arguments
 * in the coefficient and in the trigonometric part of the term are fewer or as many as the series'
 * ones.
 * Otherwise an assertion fails and the program aborts. base_pseries::merge_args,
 * base_pseries::append_cf_args and base_pseries::append_trig_args can be used to add the needed arguments
 * to the series.
 *
 * This function performs some checks and then calls base_pseries::ll_insert.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline typename base_pseries<Cf, Trig, Term, I, Derived>::it_s_index base_pseries<Cf, Trig, Term, I, Derived>::insert(
    const term_type &term, bool sign, const it_s_index *it_hint)
  {
    const size_t cw=cf_width(), tw=trig_width();
// It should not happen because resizing in this case should already be managed
// by addition and multiplication routines.
    p_assert(!term.g_cf().larger(cw));
    p_assert(!term.g_trig().larger(tw));
    term_type *new_term=0;
    if (term.g_cf().smaller(cw) || term.g_trig().smaller(tw))
    {
      new_term = new term_type(term);
      new_term->increase_size(cw,tw);
    }
    if (term.g_trig().sign()<0)
    {
      if (new_term==0)
      {
        new_term = new term_type(term);
      }
      new_term->invert_trig_args();
    }
    const term_type *insert_term;
    if (new_term==0)
    {
      insert_term=&term;
    }
    else
    {
      insert_term=new_term;
    }
    // We must check for ignorability here, since we need term's cf to be the right size when dealing
    // with non-numerical cf (i.e., polynomials, etc.) by passing the vector of symbols.
    if (insert_term->is_ignorable(cf_s_vec_))
    {
      return s_index().end();
    }
    it_s_index ret_it=ll_insert(*insert_term,sign,it_hint);
    delete new_term;
    return ret_it;
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2>
    inline typename base_pseries<Cf, Trig, Term, I, Derived>::it_s_index base_pseries<Cf, Trig, Term, I, Derived>::insert(
    const Term<Cf2,trig_type> &term, bool sign, const it_s_index *it_hint)
  {
//        BOOST_STATIC_ASSERT(sizeof(U)==0);
    return insert(term_type(term),sign,it_hint);
  }

// --------------

/// Merge arguments.
/** Merge arguments with those of ps2. If the operation succeeds the size of *this will
 * be equal or greater than ps2's. We need the template U because we want to be able to merge
 * symbols also with complex counterparts.
 * @param[in] ps2 piranha::base_pseries arguments are to be merged with.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2, class Derived2>
    inline bool base_pseries<Cf, Trig, Term, I, Derived>::merge_args(const base_pseries<Cf2,trig_type,Term,I,Derived2> &ps2)
  {
    if ((void *)this==(void *)&ps2)
    {
      std::cout << "Trying to merge with self, returning true." << std::endl;
      return true;
    }
    if (!args_compatible(ps2))
    {
      std::cout << "The base_pseries are not args_compatible." << std::endl;
      if (args_different(ps2))
      {
        std::cout << "But they are args_different. Lolrus!" << std::endl;
        prepend_cf_args(ps2.cf_s_vec());
        prepend_trig_args(ps2.trig_s_vec());
        return true;
      }
      else
      {
        std::exit(1);
        return false;
      }
    }
    size_t w1=cf_width(), w2=ps2.cf_width();
    if (w2>w1)
    {
      append_cf_args(vector_psym_p(ps2.cf_s_vec().begin()+w1,ps2.cf_s_vec().end()));
    }
    w1=trig_width();
    w2=ps2.trig_width();
    if (w2>w1)
    {
      append_trig_args(vector_psym_p(ps2.trig_s_vec().begin()+w1,ps2.trig_s_vec().end()));
    }
    return true;
  }

/// Cumulative crop.
/**
 * Crop the series from the bottom so that the sum of the norms of the cropped terms is not
 * greater than delta.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::cumulative_crop(const double &delta)
  {
// Won't crop a nil series
    if (length()==0)
    {
      return;
    }
    double part_norm=0.;
    it_s_index it=boost::prior(s_index().end());
    while (1)
    {
      part_norm+=it->norm(cf_s_vec_);
      if (part_norm>=delta)
      {
        break;
      }
      if (it==s_index().begin())
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
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::crop(const double &delta)
  {
// Won't crop a nil series
    if (length()==0)
    {
      return;
    }
    it_s_index it=boost::prior(s_index().end());
    while (1)
    {
      if (it->norm(cf_s_vec_)>=delta)
      {
        break;
      }
      if (it==s_index().begin())
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
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::crop(const it_s_index &it_f)
  {
// Won't crop a nil series or if the crop limit is the end of the series
    if (length()==0 || it_f==s_index().end())
    {
      return;
    }
    it_s_index it=boost::prior(s_index().end());
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
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::spectral_cutoff(const double &achieved_tdp,
    const double &desired_sdp)
  {
    crop(sdp_cutoff(achieved_tdp,desired_sdp));
  }
}
#endif
