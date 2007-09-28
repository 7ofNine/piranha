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

#ifndef PIRANHA_LIGHT_TERM_H
#define PIRANHA_LIGHT_TERM_H

namespace piranha
{
/// Lightweight term to be used in series multiplications.
  template <class Cf, class Trig>
    struct light_term
  {
// Manipulation.
    void invert_trig_args()
    {
      trig.invert_sign();
      if (!trig.g_flavour())
      {
// FIXME: maybe here a invert_sign function for the coefficient should be used?
        cf*=-1;
      }
    }
    bool operator==(const light_term &t) const
    {
      return (trig == t.trig);
    }
// Data members.
    mutable Cf    cf;
    mutable Trig  trig;
  };

/// Overload for the computation of the hash_value.
  template <class Cf, class Trig>
    inline size_t hash_value(const light_term<Cf,Trig> &t)
  {
    return t.trig.hasher();
  }
}

#endif
