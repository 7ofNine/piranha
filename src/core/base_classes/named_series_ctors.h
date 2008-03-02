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

#ifndef PIRANHA_NAMED_SERIES_CTORS
#define PIRANHA_NAMED_SERIES_CTORS

namespace piranha
{
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline named_series<__PIRANHA_NAMED_SERIES_TP>::named_series(const std::string &filename)
  {
    read_from_file(filename);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline named_series<__PIRANHA_NAMED_SERIES_TP>::named_series()
  {}
}

#endif
