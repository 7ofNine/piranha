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

#ifndef PIRANHA_BASE_PSERIES_HOOKS_H
#define PIRANHA_BASE_PSERIES_HOOKS_H

namespace piranha
{
/// Hooks for piranha::base_pseries.
  template <class Derived>
    class base_pseries_hooks
  {
    protected:
// Hooks.
/// Default implementation of assignment hook.
/**
 * Templatized this way because we want to be able to assign real series to complex ones.
 */
      template <class Derived2>
        void assignment_hook(const Derived2 &) {}
/// Default implementation of swap hook.
      void swap_hook(Derived &) {}
/// Default implementation of the hook for post-insertion of a new term.
      template <class Term>
        void new_term_post_insertion_hook(const Term &) {}
/// Default implementation of the hook for post-erase of a term.
      template <class Term>
        void term_pre_erase_hook(const Term &) {}
/// Default implementation of the hook for pre-update of a term.
      template <class Term, class Cf>
        void term_pre_update_hook(const Term &, const Cf &) {}
  };
}
#endif
