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

#ifndef PIRANHA_NORM_TOOLBOX_H
#define PIRANHA_NORM_TOOLBOX_H

namespace piranha
{
/// Norm toolbox.
/**
 * The norm should kept up-to-date during term insertions. It is calculated using a norm() method
 * provided by the coefficient class.
 */
  template <class Derived>
    class norm_toolbox
  {
    public:
/// Default ctor.
      norm_toolbox():norm_(0.)
        {}
/// Copy ctor.
      norm_toolbox(const norm_toolbox &n):norm_(n.norm_)
        {}
      void upgrade_norm(const double &new_real)
      {
        norm_+=std::abs(new_real);
      }
      void downgrade_norm(const double &new_real)
      {
        norm_-=std::abs(new_real);
      }
/// Const getter for norm.
      const double &g_norm() const
      {
        return norm_;
      }
    protected:
/// Getter for norm.
      double &norm()
      {
        return norm_;
      }
    protected:
      double    norm_;
  };
}

#endif
