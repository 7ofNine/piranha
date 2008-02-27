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

#ifndef PIRANHA_NTUPLE_H
#define PIRANHA_NTUPLE_H

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>

namespace piranha
{
  /// Wrapper for tuple of homogeneous types.
  template <class T, int N>
    struct ntuple
  {
    BOOST_STATIC_ASSERT(N > 0);
    typedef boost::tuples::cons<T,typename ntuple<T,N-1>::type> type;
  };
  
  template <class T>
    struct ntuple<T,1>
  {
    typedef boost::tuples::cons<T,boost::tuples::null_type> type;
  };
}

#endif
