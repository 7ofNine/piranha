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

#ifndef PIRANHA_IMONOMIAL_H
#define PIRANHA_IMONOMIAL_H

#include <boost/static_assert.hpp>

namespace piranha
{
// Forward declaration.
  template <class Cf, class Index, class Expo, class Derived>
    class base_ipoly;

/// Indexed monomial.
  template <class Cf, class Index>
    class imonomial
  {
      BOOST_STATIC_ASSERT(boost::is_integral<Index>::value);
      BOOST_STATIC_ASSERT(boost::is_pod<Index>::value);
      BOOST_STATIC_ASSERT(!(boost::integer_traits<Index>::is_signed));
      template <class Cf2, class Index2, class Expo, class Derived>
        friend class base_ipoly;
    public:
      imonomial(const Cf &value, const Index &i):index(i),cf(value)
        {}
      ~imonomial()
        {}
      bool operator==(const imonomial &m) const
      {
        return (index == m.index);
      }
// TODO: needed?
      bool operator<(const imonomial &m) const
      {
        return (index < m.index);
      }
      const Cf &g_cf() const
      {
        return cf;
      }
      const Index &g_index() const
      {
        return index;
      }
    private:
      imonomial()
        {}
    private:
      Index               index;
      mutable Cf          cf;
    public:
      static const Index  max_index = ((((((Index)1)<<(sizeof(Index)*8-1))-1)<<1)+1);
  };

  template <class Cf, class Index>
    const Index imonomial<Cf,Index>::max_index;

// Hasher functor to be used in polynomial container.
  template <class Monomial>
    struct monomial_hasher
  {
    const size_t &operator()(const Monomial &m) const
    {
      return m.g_index();
    }
  };
}

#endif
