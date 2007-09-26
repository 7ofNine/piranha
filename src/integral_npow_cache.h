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

#ifndef PIRANHA_INTEGRAL_NPOW_CACHE_H
#define PIRANHA_INTEGRAL_NPOW_CACHE_H

#include <boost/integer_traits.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_pod.hpp>
#include <cmath>
#include <vector>

#include "p_assert.h"

namespace piranha
{
/// Natural power cacher for integral type.
  template <class T>
    class integral_npow_cache
  {
      BOOST_STATIC_ASSERT(boost::is_integral<T>::value);
      BOOST_STATIC_ASSERT(boost::is_pod<T>::value);
      BOOST_STATIC_ASSERT(!(boost::integer_traits<T>::is_signed));
      typedef std::vector<std::vector<T> > container;
    public:
      static const T &request(const int &n, const T &arg)
      {
        p_assert(n >= 0);
        p_assert(arg >= 0);
        while (arg >= private_cache_.size())
        {
// Add a row to the matrix.
          private_cache_.push_back(std::vector<T>());
// Add the first element to the row.
          private_cache_.back().push_back(1);
        }
        while ((size_t)n >= private_cache_[arg].size())
        {
//std::cout << "before: " << private_cache_[arg].back() << '\n';
          private_cache_[arg].push_back(arg*private_cache_[arg].back());
//std::cout << "now: " << private_cache_[arg].back() << '\n';
        }
        return private_cache_[arg][n];
      }
    private:
      static container    private_cache_;
  };

  template <class T>
    typename integral_npow_cache<T>::container integral_npow_cache<T>::private_cache_;
}

#endif
