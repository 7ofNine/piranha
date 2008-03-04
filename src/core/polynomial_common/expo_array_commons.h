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

#ifndef PIRANHA_EXPO_ARRAY_COMMONS_H
#define PIRANHA_EXPO_ARRAY_COMMONS_H

#include <boost/static_assert.hpp>
#include <iostream>

#define derived_const_cast (static_cast<Derived const *>(this))
#define derived_cast (static_cast<Derived *>(this))

namespace piranha
{
  /// Common class for dense exponent array.
  /**
   * Intended to add specific methods to plain arrays for the manipulation of exponent
   * parts in polynomials.
   */
  template <class Derived, int Pos>
    class expo_array_commons
  {
      BOOST_STATIC_ASSERT(Pos >= 0);
    public:
      // I/O.
      template <class ArgsTuple>
        void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const
      {
        // We assert like this because we want to make sure we don't go out of boundaries,
        // and because in case of fixed-width we may have smaller size of v wrt to "real" size.
        p_assert(args_tuple.template get<Pos>().size() <= derived_const_cast->size());
        derived_const_cast->print(out_stream);
      }
      void print_latex(std::ostream &out_stream, const vector_psym_p &v) const
      {
        // TODO: implement.
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
