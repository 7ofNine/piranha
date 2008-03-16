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

#ifndef PIRANHA_CF_SERIES_H
#define PIRANHA_CF_SERIES_H

#include <iostream>
#include <string>

#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_CF_SERIES_TP_DECL class Derived
#define __PIRANHA_CF_SERIES_TP Derived

namespace piranha
{
  /// Toolbox for using a series as a coefficient in another series.
  /**
   * Intended to be inherited by piranha::base_series.
   */
  template <__PIRANHA_CF_SERIES_TP_DECL>
    class cf_series
  {
    public:
      template <class ArgsTuple>
        void print_plain(std::ostream &, const ArgsTuple &) const;
      template <class ArgsTuple>
        bool is_insertable(const ArgsTuple &) const;
      template <class ArgsTuple>
        bool needs_padding(const ArgsTuple &) const;
      template <class ArgsTuple>
        bool is_ignorable(const ArgsTuple &) const;
      template <class ArgsTuple>
        void pad_right(const ArgsTuple &);
      template <class ArgsTuple>
        Derived &add(const Derived &, const ArgsTuple &);
      template <class ArgsTuple>
        Derived &subtract(const Derived &, const ArgsTuple &);
      template <class ArgsTuple>
        void invert_sign(const ArgsTuple &);
      void swap(Derived &);
      template <class ArgsTuple, class Layout>
        void apply_layout(const ArgsTuple &, const Layout &);
      template <class ArgsTuple>
        Derived &mult_by(const int &, const ArgsTuple &);
      template <class ArgsTuple>
        Derived &mult_by(const double &, const ArgsTuple &);
    protected:
      template <class ArgsTuple>
        void construct_from_string(const std::string &, const ArgsTuple &);
      template <class T, class ArgsTuple>
        Derived &mult_by_generic(const T &, const ArgsTuple &);
  };
}

#include "cf_series_io.h"
#include "cf_series_manip.h"
#include "cf_series_math.h"
#include "cf_series_probe.h"

#undef __PIRANHA_CF_SERIES_TP
#undef __PIRANHA_CF_SERIES_TP_DECL
#undef derived_const_cast
#undef derived_cast

#endif
