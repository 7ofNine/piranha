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

#include <boost/static_assert.hpp>
#include <memory>

#include "../p_assert.h"
#include "../utils.h"                             // For class_converter.

// Useful shortcuts.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_BASE_SERIES_TP_DECL class Term, class Derived, class Allocator
#define __PIRANHA_BASE_SERIES_TP Term,Derived,Allocator

namespace piranha
{
  /// Base series class.
  /**
   * Term must derive from piranha::base_term class.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL> class base_series
  {
      /// Alias for term type.
      typedef Term term_type;
      /// Alias for coefficient type.
      typedef typename term_type::cf_type cf_type;
      /// Alias for key type.
      typedef typename term_type::key_type key_type;
      /// Alias for allocator type.
      typedef Allocator allocator_type;
      /// Alias for allocator rebinding to term_type.
      typedef typename allocator_type::template rebind<term_type>::other term_allocator_type;
      BOOST_STATIC_ASSERT(term_type::size == 2);
    public:
      template <bool, bool, class Term2, class ArgsTuple, class SortedIterator>
        SortedIterator insert(const Term2 &, const ArgsTuple &, SortedIterator);
      template <bool, bool, class Term2, class SortedIterator>
        SortedIterator insert(const Term2 &, SortedIterator);
      template <class Term2, class ArgsTuple, class SortedIterator>
        SortedIterator insert(const Term2 &, const ArgsTuple &, SortedIterator);
      template <class Term2, class SortedIterator>
        SortedIterator insert(const Term2 &, SortedIterator);
      template <int N, class ArgsTuple, class Iterator>
        void term_erase(const ArgsTuple &, Iterator);
      template <int N, class Iterator>
        void term_erase(Iterator);
    private:
      template <class PinpointIterator>
        PinpointIterator find_term(const term_type &) const;
      template <bool, class ArgsTuple, class SortedIterator>
        SortedIterator ll_insert(const term_type &, const ArgsTuple &, SortedIterator);
      template <bool, class ArgsTuple, class SortedIterator>
        SortedIterator term_insert_new(const term_type &, const ArgsTuple &, SortedIterator);
      template <class ArgsTuple, class PinpointIterator>
        void term_update(const ArgsTuple &, PinpointIterator, cf_type &);
      // Functors.
      struct modifier_invert_term_sign
      {
        void operator()(term_type &term) const
        {
          term.elements.template get<0>().invert_sign();
        }
      };
      struct modifier_update_cf
      {
          modifier_update_cf(cf_type &new_cf):m_new_cf(new_cf) {}
          ~modifier_update_cf() {}
          void operator()(term_type &term)
          {
            term.cf().swap(m_new_cf);
          }
        private:
          cf_type &m_new_cf;
      };
    protected:
      /// Rebound allocator for term type.
      static term_allocator_type term_allocator;
  };

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    typename base_series<__PIRANHA_BASE_SERIES_TP>::term_allocator_type
    base_series<__PIRANHA_BASE_SERIES_TP>::term_allocator;
}

#include "base_series_manip.h"

#undef derived_const_cast
#undef derived_cast

#endif
