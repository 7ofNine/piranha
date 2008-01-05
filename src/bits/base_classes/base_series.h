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
#include "../utils.h" // For class_converter.
#define derived_cast static_cast<Derived const *>(this)

namespace piranha
{
  // Base series class.
  /*
   * Term must derive from piranha::base_term class. 
   */
  template <class Term, class Derived> class base_series
  {
    public:
      template <class SortedIterator, class Term2, bool AdditionalChecks, bool Sign> SortedIterator insert(
        const Term2 &term_, SortedIterator it_hint)
      {
        typedef typename Derived::term_type term_type;
        class_converter<Term2,term_type> term(term_);
        // It should not happen because resizing in this case should already be managed
        // by external routines (merge_args, and input from file).
        p_assert(term.result.is_insertable(*this));
        term_type *new_term(0);
        const bool padding_needed=(term.result.needs_padding(*this));
        if (unlikely(padding_needed))
        {
          new_term=derived_cast->term_allocator.allocate(1);
          term_allocator.construct(new_term, term.result);
          new_term->pad_right(*this);
        }
        if (AdditionalChecks)
        {
          derived_cast->i_perform_additional_checks(term.result,new_term);
        }
        const term_type *insert_term;
        if (new_term == 0)
        {
          insert_term=&term.result;
        } else
        {
          insert_term=new_term;
        }
        SortedIterator ret_it=derived_cast->ll_insert<Sign>(*insert_term, it_hint);
        if (new_term != 0)
        {
          derived_cast->term_allocator.destroy(new_term);
          derived_cast->term_allocator.deallocate(new_term, 1);
        }
        return ret_it;
      }
  };
}

#undef derived_cast

#endif
