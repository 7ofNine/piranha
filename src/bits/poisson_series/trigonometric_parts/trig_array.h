/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_TRIG_ARRAY_H
#define PIRANHA_TRIG_ARRAY_H

#include <boost/integer_traits.hpp>
#include <memory>

#include "../base_classes/trig_array_commons.h"
#include "../../base_classes/int_array.h"

namespace piranha
{
  template <int Bits, class Allocator = std::allocator<char> >
/// Trigonometric array, dynamically sized version.
/**
 * It wraps a piranha::int_array with signed integer sized Bits, and adds the
 * capabilities needed for trigonometric manipulation.
 */
    class trig_array: public int_array<Bits,true,Allocator>,
      public trig_array_commons<trig_array<Bits> >
  {
      typedef trig_array_commons<trig_array> trig_commons;
      typedef int_array<Bits,true,Allocator> ancestor;
    public:
      typedef typename ancestor::value_type value_type;
      typedef typename ancestor::size_type size_type;
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_array():ancestor::int_array() {}
/// Ctor from piranha::deque_string.
      trig_array(const deque_string &sd):ancestor::int_array(),trig_commons::trig_array_commons(sd) {}
/// Pad right.
      void pad_right(const size_t &n)
      {
        p_assert(n >= ancestor::size());
        ancestor::resize(n);
      }
// Probing.
/// Data footprint.
/**
 * Returns the memory occupied by the data members.
 */
      size_t data_footprint() const {return (ancestor::size()*sizeof(value_type));}
      template <class Series>
        bool checkup(const Series &s) const
      {
        switch (s.trig_width() != ancestor::size())
        {
          case true:
            std::cout << "Size mismatch in trig_array." << std::endl;
            return false;
          default:
            return true;
        }
      }
/// Does trig_array needs padding to be inserted in series with trig_width equal to n?
      bool needs_padding(const size_t &n) const {return (ancestor::size() < n);}
/// Does is trig_array insertable in series with trig_width equal to n?
      bool is_insertable(const size_t &n) const {return (ancestor::size() <= n);}
// End INTERFACE definition.
//-------------------------------------------------------
  };
}
#endif
