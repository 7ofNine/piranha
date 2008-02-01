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

#ifndef PIRANHA_BASE_SERIES_H
#define PIRANHA_BASE_SERIES_H

#include "../p_assert.h"
#include "../utils.h"                             // For class_converter.

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  /// Base series class.
  /**
   * Term must derive from piranha::base_term class.
   */
  template <class Term, class Derived> class base_series
  {
    public:
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
      template <class SortedIterator, class Term2, bool AdditionalChecks, bool Sign> SortedIterator insert(
        const Term2 &term_, SortedIterator it_hint)
      {
        class_converter<Term2,Term> term(term_);
        // It should not happen because resizing in this case should already be managed
        // by external routines (merge_args, and input from file).
        p_assert(term.result.is_insertable(*derived_const_cast));
        Term *new_term(0);
        const bool padding_needed=(term.result.needs_padding(*derived_const_cast));
        switch(unlikely(padding_needed))
        {
          case true:
            new_term=Derived::term_allocator.allocate(1);
            Derived::term_allocator.construct(new_term, term.result);
            new_term->pad_right(*derived_const_cast);
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
        SortedIterator ret_it=derived_cast->template ll_insert<Sign>(*insert_term, &it_hint);
        switch (new_term == 0)
        {
          case true:
            break;
          case false:
            Derived::term_allocator.destroy(new_term);
            Derived::term_allocator.deallocate(new_term, 1);
        }
        return ret_it;
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
