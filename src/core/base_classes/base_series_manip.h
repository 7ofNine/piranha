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

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp> // For iterator type detection.

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
    template <bool CanonicalCheck, bool Sign, class Term2, class ArgsTuple, class SortedIterator>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term_,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    // We need to do this because when doing insert_new we may need to change sign. We need a non-const sorted
    // iterator to do that.
    BOOST_STATIC_ASSERT((boost::is_same<SortedIterator,typename Derived::sorted_iterator>::value));
    class_converter<Term2,term_type> converted_term(term_);
    // Make sure the appropriate routines for the management of arguments have been called.
    p_assert(converted_term.result.is_insertable(args_tuple));
    term_type *new_term(0);
    switch (unlikely(converted_term.result.needs_padding(args_tuple)))
    {
      case true:
        new_term=term_type::allocator.allocate(1);
        term_type::allocator.construct(new_term,converted_term.result);
        new_term->pad_right(args_tuple);
        break;
      case false:
        ;
    }
    if (CanonicalCheck)
    {
      if (!converted_term.result.is_canonical(args_tuple))
      {
        if (new_term == 0)
        {
          new_term=term_type::allocator.allocate(1);
          term_type::allocator.construct(new_term,converted_term.result);
        }
        new_term->canonicalise(args_tuple);
      }
    }
    const term_type *insert_term(0);
    switch (new_term == 0)
    {
      case true:
        insert_term=&converted_term.result;
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

  // This cannot be const because we use this in insertion function, hence we need a non const iterator.
  /// Find term.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class PinpointIterator>
    inline PinpointIterator base_series<__PIRANHA_BASE_SERIES_TP>::find_term(const term_type &t)
  {
    BOOST_STATIC_ASSERT((boost::is_same<PinpointIterator,typename Derived::pinpoint_iterator>::value));
    return derived_cast->template nth_index<1>().find(t);
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <bool Sign, class ArgsTuple, class SortedIterator>
    inline SortedIterator base_series<__PIRANHA_BASE_SERIES_TP>::ll_insert(const term_type &term,
    const ArgsTuple &args_tuple, SortedIterator it_hint)
  {
    BOOST_STATIC_ASSERT((boost::is_same<SortedIterator,typename Derived::sorted_iterator>::value));
    typedef typename Derived::pinpoint_iterator pinpoint_iterator;
    if (term.is_ignorable(args_tuple))
    {
      return derived_const_cast->template nth_index<0>().end();
    }
    p_assert(term.is_insertable(args_tuple) and !term.needs_padding(args_tuple) and term.is_canonical(args_tuple));
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
          new_c.add(term.m_cf,args_tuple);
          break;
        case false:
          new_c=it->m_cf;
          new_c.subtract(term.m_cf,args_tuple);
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
    BOOST_STATIC_ASSERT((boost::is_same<SortedIterator,typename Derived::sorted_iterator>::value));
    typename arg_manager<Term>::arg_assigner aa(args_tuple);
    SortedIterator it_new(derived_cast->template nth_index<0>().insert(it_hint,term));
    // TODO: use asserts here? The problem here is that we are using hinted
    // insertion, the return value is different from above (but above an assert
    // is needed too).
    p_assert(it_new != derived_const_cast->template nth_index<0>().end());
    if (!Sign)
    {
      // This is an O(1) operation, since the order in the set is not changed
      // There is a re-hash involved, it still should be cheaper than
      // creating a new term though.
      modifier_invert_term_sign<ArgsTuple> m(args_tuple);
      action_assert(derived_cast->template nth_index<0>().modify(it_new,m));
    }
    return it_new;
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <int N, class ArgsTuple, class Iterator>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::term_erase(const ArgsTuple &args_tuple,
    Iterator it)
  {
    typename arg_manager<Term>::arg_assigner aa(args_tuple);
    derived_cast->template nth_index<N>().erase(it);
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple, class PinpointIterator>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::term_update(const ArgsTuple &args_tuple,
    PinpointIterator it, cf_type &new_c)
  {
    BOOST_STATIC_ASSERT((boost::is_same<PinpointIterator,typename Derived::pinpoint_iterator>::value));
    typename arg_manager<Term>::arg_assigner aa(args_tuple);
    // Update the existing term.
    modifier_update_cf m(new_c);
    action_assert(derived_cast->template nth_index<1>().modify(it,m));
  }

  /// Swap the terms with another series.
  /**
   * All terms get swapped.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::swap_terms(Derived &ps2)
  {
    p_assert(derived_cast != &ps2);
    derived_cast->m_container.swap(ps2.m_container);
  }

  /// Apply an arguments layout to all terms and insert them into retval.
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple, class Layout>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::apply_layout_to_terms(
    const ArgsTuple &args_tuple, const Layout &l, Derived &retval) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    typedef typename Derived::sorted_iterator sorted_iterator;
    typedef typename Derived::term_type term_type;
    const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
    sorted_iterator it_hint = retval.template nth_index<0>().end();
    for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin();
      it != it_f; ++it)
    {
      term_type term(*it);
      term.m_cf.apply_layout(args_tuple,l);
      term.m_key.apply_layout(args_tuple,l);
      it_hint = retval.insert(term,args_tuple,it_hint);
    }
  }
}

#endif