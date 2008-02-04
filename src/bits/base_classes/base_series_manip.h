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

#ifndef PIRANHA_BASE_SERIES_MANIP_H
#define PIRANHA_BASE_SERIES_MANIP_H

#include "../stats.h"

namespace piranha
{
  // TODO: update doc here.
  /// High-level insertion function.
  /**
   * This function is used to insert terms into a series. It requires that the number of arguments
   * of each element of the term is smaller or equal to the series',
   * otherwise an assertion fails and the program aborts. base_pseries::merge_args,
   * base_pseries::append_cf_args, base_pseries::append_trig_args, etc. can be used to add the needed arguments
   * to the series.
   *
   * This function performs some checks and then calls base_pseries::ll_insert.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class Term2, class ArgsTuple, class SortedIterator, bool AdditionalChecks, bool Sign>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term_,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    class_converter<Term2,Term> term(term_);
    // It should not happen because resizing in this case should already be managed
    // by other routines (merge_args, input from file, etc.).
    p_assert(term.result.is_insertable(args_tuple));
    Term *new_term(0);
    switch(unlikely(term.result.needs_padding(args_tuple)))
    {
      case true:
        new_term=term_allocator.allocate(1);
        term_allocator.construct(new_term, term.result);
        new_term->pad_right(args_tuple);
        break;
      case false:
        ;
    }
    if (AdditionalChecks)
    {
      new_term = derived_const_cast->i_perform_additional_checks(term.result,new_term);
    }
    const Term *insert_term(0);
    switch (new_term == 0)
    {
      case true:
        insert_term=&term.result;
        break;
      case false:
        insert_term=new_term;
    }
    SortedIterator ret_it=ll_insert<Sign>(*insert_term,args_tuple,it_hint);
    switch (new_term == 0)
    {
      case true:
        break;
      case false:
        term_allocator.destroy(new_term);
        term_allocator.deallocate(new_term,1);
    }
    return ret_it;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class Term2, class SortedIterator, bool AdditionalChecks, bool Sign>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term,
    SortedIterator it_hint)
  {
    return insert<Term2,typename Derived::arguments_tuple_type,SortedIterator,AdditionalChecks,Sign>(
      term,derived_const_cast->m_arguments,it_hint);
  }

  /// Find term.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class PinpointIterator>
    inline PinpointIterator base_series<__PIRANHA_BASE_SERIES_TP>::find_term(const term_type &t) const
  {
    return derived_const_cast->g_p_index().find(t);
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool Sign, class ArgsTuple, class SortedIterator>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::ll_insert(const term_type &term,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    typedef typename Derived::pinpoint_iterator pinpoint_iterator;
    // We must check for ignorability here, before assertions, since at this point cf_width could be zero
    // for polynomials, and thus assertion would wrongly fail.
    if (term.is_ignorable(args_tuple))
    {
      return derived_const_cast->g_s_index().end();
    }
    p_assert(term.is_insertable(args_tuple) and !term.needs_padding(args_tuple));
    //p_assert(term.trig().sign()>0);
    SortedIterator ret_it;
    pinpoint_iterator it(find_term<pinpoint_iterator>(term));
    if (it == derived_const_cast->g_p_index().end())
    {
      // The term is NOT a duplicate, insert in the set. Record where we inserted,
      // so it can be used in additions and multiplications.
      ret_it=derived_cast->template term_insert_new<Sign>(term,it_hint);
      stats::insert();
    }
    else
    {
      // The term is in the set, hence an existing term will be modified.
      // Add or subtract according to request.
      cf_type new_c;
      switch (Sign)
      {
        case true:
          new_c=it->elements.template get<0>();
          new_c.add(term.elements.template get<0>());
          break;
        case false:
          new_c=it->elements.template get<0>();
          new_c.subtract(term.elements.template get<0>());
      }
      // Check if the resulting coefficient can be ignored (ie it is small).
      if (new_c.is_ignorable(args_tuple))
      {
        derived_cast->term_erase(it);
      }
      else
      {
        derived_cast->term_update(it,new_c);
      }
      // If we are erasing or updating there's no point in giving an hint on where
      // the action took place, just return the end() iterator.
      ret_it=derived_cast->g_s_index().end();
      stats::pack();
    }
    return ret_it;
  }
}

#endif
