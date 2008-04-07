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
#include "../type_traits.h"
#include "../utils.h" // For class_converter.

// Useful shortcuts.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_BASE_SERIES_TP_DECL class Term, char Separator, class Allocator, class Derived
#define __PIRANHA_BASE_SERIES_TP Term,Separator,Allocator,Derived

namespace piranha
{
  /// Base series class.
  /**
   * Term must derive from piranha::base_term class.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    class base_series
  {
      /// Alias for term type.
      typedef Term term_type;
      /// Alias for coefficient type.
      typedef typename term_type::cf_type cf_type;
      /// Alias for key type.
      typedef typename term_type::key_type key_type;
      /// Alias for allocator type.
      typedef Allocator allocator_type;
      // Evaluation type. Used internally.
      typedef typename eval_type<Derived>::type eval_type;
    public:
      template <bool, bool, class Term2, class ArgsTuple, class SortedIterator>
        SortedIterator insert(const Term2 &, const ArgsTuple &, SortedIterator);
      template <class Term2, class ArgsTuple, class SortedIterator>
        SortedIterator insert(const Term2 &, const ArgsTuple &, SortedIterator);
      template <class ArgsTuple>
        double calculate_norm(const ArgsTuple &) const;
    protected:
      static const char separator = Separator;
      // Check that the separators do not conflict.
      BOOST_STATIC_ASSERT(separator != term_type::separator);
      template <class Number, class ArgsTuple>
        void construct_from_number(const Number &, const ArgsTuple &);
      template <class ArgsTuple>
        void print_terms_plain(std::ostream &, const ArgsTuple &, int limit) const;
      template <class ArgsTuple>
        void print_terms_latex(std::ostream &, const ArgsTuple &, int limit) const;
      template <int N, class ArgsTuple, class Iterator>
        void term_erase(const ArgsTuple &, Iterator);
      void swap_terms(Derived &);
      template <class ArgsTuple, class Layout>
        void apply_layout_to_terms(const ArgsTuple &, const Layout &, Derived &) const;
      template <bool, class Derived2, class ArgsTuple>
        void merge_terms(const Derived2 &, const ArgsTuple &);
      template <class T, class ArgsTuple>
        void multiply_coefficients_by(const T &, Derived &, const ArgsTuple &) const;
      template <class T, class ArgsTuple>
        void divide_coefficients_by(const T &, Derived &, const ArgsTuple &) const;
      template <bool, class Number, class ArgsTuple>
        Derived &merge_with_number(const Number &, const ArgsTuple &);
      bool is_single_cf() const;
      template <class ArgsTuple>
        eval_type time_evaluation(const double &, const ArgsTuple &) const;
    private:
      template <class PinpointIterator>
        PinpointIterator find_term(const term_type &);
      template <bool, class ArgsTuple, class SortedIterator>
        SortedIterator ll_insert(const term_type &, const ArgsTuple &, SortedIterator);
      template <bool, class ArgsTuple, class SortedIterator>
        SortedIterator term_insert_new(const term_type &, const ArgsTuple &, SortedIterator);
      template <class ArgsTuple, class PinpointIterator>
        void term_update(const ArgsTuple &, PinpointIterator, cf_type &);
      // Functors.
      template <class ArgsTuple>
        struct modifier_invert_term_sign
      {
        modifier_invert_term_sign(const ArgsTuple &args_tuple):a(args_tuple) {}
        // This (and below) are passed as const to allow modification from const interators
        // that are often used in standard implementations of data structures.
        void operator()(const term_type &term) const
        {
          term.m_cf.invert_sign(a);
        }
        const ArgsTuple &a;
      };
      struct modifier_update_cf
      {
          modifier_update_cf(cf_type &new_cf):m_new_cf(new_cf) {}
          ~modifier_update_cf() {}
          void operator()(const term_type &term)
          {
            term.m_cf.swap(m_new_cf);
          }
        private:
          cf_type &m_new_cf;
      };
  };
}

#include "base_series_io.h"
#include "base_series_manip.h"
#include "base_series_math.h"
#include "base_series_probe.h"

#undef derived_const_cast
#undef derived_cast
#undef __PIRANHA_BASE_SERIES_TP_DECL
#undef __PIRANHA_BASE_SERIES_TP

#endif
