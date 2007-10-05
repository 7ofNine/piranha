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

#ifndef PIRANHA_HASH_SET_HM_H
#define PIRANHA_HASH_SET_HM_H

#include <ext/hash_set>

namespace piranha
{
  template <class Element, class Hasher, class Eq, class Allocator, bool StoreHash = false>
    class mult_hash
  {
      typedef __gnu_cxx::hash_set<Element,Hasher,Eq,Allocator> container_type;
    public:
      typedef typename container_type::iterator iterator;
      mult_hash()
        {}
      mult_hash(const size_t &s):private_container_(s)
        {}
      iterator begin() const
      {
        return private_container_.begin();
      }
      iterator end() const
      {
        return private_container_.end();
      }
      iterator find(const Element &e) const
      {
        return private_container_.find(e);
      }
      iterator insert(const Element &e)
      {
        return private_container_.insert(e).first;
      }
    private:
      __gnu_cxx::hash_set<Element,Hasher,Eq,Allocator>    private_container_;
  };
}

#endif
