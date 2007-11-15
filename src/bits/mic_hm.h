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

#ifndef PIRANHA_MIC_HM_H
#define PIRANHA_MIC_HM_H

#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index_container.hpp>

namespace piranha
{
  template <class Element, class Hasher, class Eq, class Allocator, bool StoreHash>
    class mult_hash
  {
      typedef boost::multi_index_container<
        Element,
        boost::multi_index::indexed_by<
          boost::multi_index::hashed_unique<boost::multi_index::identity<Element> >
        >,
      Allocator> container_type;
    public:
      typedef typename container_type::const_iterator iterator;
      typedef typename container_type::iterator point_iterator;
      mult_hash()
        {
          private_container_.max_load_factor(.8);
        }
      mult_hash(const size_t &)
        {
          private_container_.max_load_factor(.8);
        }
      iterator begin() const
      {
        return private_container_.begin();
      }
      iterator end() const
      {
        return private_container_.end();
      }
      point_iterator find(const Element &e) const
      {
        return private_container_.find(e);
      }
      point_iterator insert(const Element &e)
      {
        return private_container_.insert(e).first;
      }
    private:
      container_type    private_container_;
  };

  template <class Element, class Hasher, class Eq, class Allocator, bool StoreHash>
    struct c_mult_hash
  {
    typedef mult_hash<Element,Hasher,Eq,Allocator,StoreHash> type;
  };
}

#endif
