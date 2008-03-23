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

#ifndef PIRANHA_TRUNCATORS_H
#define PIRANHA_TRUNCATORS_H

namespace piranha
{
  template <class Series1, class Series2, class ArgsTuple>
    class no_truncation
  {
      typedef typename Series1::term_type term_type;
      typedef typename Series1::term_type::cf_type cf_type1;
      typedef typename Series2::term_type::cf_type cf_type2;
      typedef typename Series1::term_type::key_type key_type;
    public:
      no_truncation(const Series1 &, const Series2 &, const ArgsTuple &) {}
      bool accept_term(const term_type &) const {return true;}
      bool skip_from_here(const cf_type1 &, const key_type &, const cf_type2 &, const key_type &) const {return false;}
  };
}

#endif
