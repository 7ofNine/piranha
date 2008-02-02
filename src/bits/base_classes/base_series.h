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

// Useful shortcuts.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_BASE_SERIES_TP_DECL class Term, class Derived
#define __PIRANHA_BASE_SERIES_TP Term,Derived

namespace piranha
{
  /// Base series class.
  /**
   * Term must derive from piranha::base_term class.
   */
  template <__PIRANHA_BASE_SERIES_TP_DECL> class base_series
  {
    public:
      template <class Term2, class ArgsTuple, class SortedIterator, bool AdditionalChecks, bool Sign>
        SortedIterator insert(const Term2 &, const ArgsTuple &, SortedIterator);
      template <class Term2, class SortedIterator, bool AdditionalChecks, bool Sign>
        SortedIterator insert(const Term2 &, SortedIterator);
  };
}

#include "base_series_manip.h"

#undef derived_const_cast
#undef derived_cast

#endif
