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

#ifndef PIRANHA_PS_TERM_MANIP_H
#define PIRANHA_PS_TERM_MANIP_H

namespace piranha
  {
  /// Invert the sign of trigonometric arguments multiplicers.
  template <class Cf,class Trig>
  inline void ps_term<Cf,Trig>::invert_trig_args()
  {
    trig_args_.invert_sign();
    if (!flavour_)
      {
        // FIXME: maybe here a invert_sign function for the coefficient should be used?
        c_*=-1;
      }
  }


  /// Resize term.
  template <class Cf,class Trig>
  inline void ps_term<Cf,Trig>::resize(const size_t &cw, const size_t &tw)
  {
    c_.resize(cw);
    trig_args_.resize(tw);
  }
}
#endif
