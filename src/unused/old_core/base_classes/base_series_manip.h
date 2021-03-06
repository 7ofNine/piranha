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
   * This function performs some checks and then calls ll_insert.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool AdditionalChecks, bool Sign, class Term2, class ArgsTuple, class SortedIterator>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term_,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    class_converter<Term2,term_type> term(term_);
    // It should not happen because resizing in this case should already be managed
    // by other routines (merge_args, input from file, etc.).
    p_assert(term.result.is_insertable(args_tuple));
    term_type *new_term(0);
    switch(unlikely(term.result.needs_padding(args_tuple)))
    {
      case true:
        new_term=term_type::allocator.allocate(1);
        term_type::allocator.construct(new_term, term.result);
        new_term->pad_right(args_tuple);
        break;
      case false:
        ;
    }
    if (AdditionalChecks)
    {
      new_term = derived_const_cast->canonicalise(term.result,new_term);
    }
    const term_type *insert_term(0);
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
        term_type::allocator.destroy(new_term);
        term_type::allocator.deallocate(new_term,1);
    }
    return ret_it;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class Term2, class ArgsTuple, class SortedIterator>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    return insert<true,true>(term,args_tuple,it_hint);
  }

  /// Find term.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class PinpointIterator>
    inline PinpointIterator base_series<__PIRANHA_BASE_SERIES_TP>::find_term(const term_type &t) const
  {
    return derived_const_cast->template nth_index<1>().find(t);
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
      return derived_const_cast->template nth_index<0>().end();
    }
    p_assert(term.is_insertable(args_tuple) and !term.needs_padding(args_tuple));
    //p_assert(term.trig().sign()>0);
    SortedIterator ret_it;
    pinpoint_iterator it(find_term<pinpoint_iterator>(term));
    if (it == derived_const_cast->template nth_index<1>().end())
    {
      // The term is NOT a duplicate, insert in the set. Record where we inserted,
      // so it can be used in additions and multiplications.
      ret_it=term_insert_new<Sign>(term,args_tuple,it_hint);
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
          new_c=it->m_cf;
          new_c.add(term.m_cf);
          break;
        case false:
          new_c=it->m_cf;
          new_c.subtract(term.m_cf);
      }
      // Check if the resulting coefficient can be ignored (ie it is small).
      if (new_c.is_ignorable(args_tuple))
      {
        term_erase<1>(args_tuple,it);
      }
      else
      {
        term_update(args_tuple,it,new_c);
      }
      // If we are erasing or updating there's no point in giving an hint on where
      // the action took place, just return the end() iterator.
      ret_it=derived_cast->template nth_index<0>().end();
      stats::pack();
    }
    return ret_it;
  }

  // Insert a new term into the series
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool Sign, class ArgsTuple, class SortedIterator>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::term_insert_new(const term_type &term,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    typename arg_manager<ArgsTuple>::arg_assigner aa(args_tuple);
    SortedIterator it_new(derived_cast->s_s_index().insert(it_hint,term));
    // TODO: use asserts here? The problem here is that we are using hinted
    // insertion, the return value is different from above (but above an assert
    // is needed too).
    p_assert(it_new!=derived_const_cast->end());
    if (!Sign)
    {
      // This is an O(1) operation, since the order in the set is not changed
      // There is a re-hash involved, it still should be cheaper than
      // creating a new term though.
      action_assert(derived_cast->s_s_index().modify(it_new,modifier_invert_term_sign()));
    }
    return it_new;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <int N, class ArgsTuple, class Iterator>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::term_erase(const ArgsTuple &args_tuple,
    Iterator it)
  {
    typename arg_manager<ArgsTuple>::arg_assigner aa(args_tuple);
    derived_cast->template nth_index<N>().erase(it);
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple, class PinpointIterator>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::term_update(const ArgsTuple &args_tuple,
    PinpointIterator it, cf_type &new_c)
  {
    typename arg_manager<ArgsTuple>::arg_assigner aa(args_tuple);
    // Update the existing term.
    action_assert(derived_cast->template nth_index<1>().modify(it,modifier_update_cf(new_c)));
  }

  /// Append an argument.
  /**
   * The argument will be appended at the end of the corresponding argument vector. This method can be used only
   * on empty series.
   * @param[in] N arguments position.
   * @param[in] p pirahna::psym_p to be appended.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <int N>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::append_arg(const psym_p p)
  {
    hard_assert(derived_const_cast->empty());
    try
    {
      // TODO: port.
      if (derived_const_cast->arg_index<N>(p->name()).first)
      {
        std::cout << "Symbol already present in series, I won't add it." << std::endl;
        return;
      }
      if (!cf_type::template admits_size<N>(derived_const_cast->arguments().template get<N>().size()+1))
      {
        throw exceptions::add_arguments(std::string("Cannot append cf argument, max_size reached.");
      }
      if (!key_type::template admits_size<N>(derived_const_cast->arguments().template get<N>().size()+1))
      {
        throw exceptions::add_arguments(std::string("Cannot append key argument, max_size reached.");
      }
    }





    typedef typename boost::tuples::element<N,typename term_type::tuple_type>::type element_type;
    try
    {
      Derived retval;
      retval.lin_args()=lin_args();
      retval.arguments()=arguments();
      const std::string description(psymbol_descr<N>());
      if (arg_index<N>(p->name()).first)
      {
        std::cout << description << " already present in series, I won't add it." << std::endl;
        return;
      }
      // Append psymbol.
      // NOTICE: maybe we can place this exception throwing also in pad_right() methods?
      if (element_type::max_size < retval.arguments().template get<N>().size()+1)
      {
        throw exceptions::add_arguments(std::string("Cannot append further ")+
          description+", max_size reached.");
      }
      retval.arguments().template get<N>().push_back(p);
      // Deal with lin_args if argument is trigonometric.
      if (N == psymbol::trig)
      {
        retval.lin_args().push_back(0);
      }
      // NOTICE: this can be probably split in 2 here if we want to use it in generic_series routines.
      const it_h_index it_f=g_h_index().end();
      it_s_index it_hint = retval.g_s_index().end();
      for (it_h_index it=g_h_index().begin();it!=it_f;++it)
      {
        // NOTICE: find a way to avoid resizes here?
        term_type tmp_term=(*it);
        tmp_term.elements.template get<N>().pad_right(retval.arguments());
        it_hint = retval.insert(tmp_term,it_hint);
      }
      swap(retval);
    }
    catch (exceptions::add_arguments &e)
      {}
  }
}

#endif
