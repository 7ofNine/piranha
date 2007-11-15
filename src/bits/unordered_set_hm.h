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

#ifndef PIRANHA_UNORDERED_SET_HM_H
#define PIRANHA_UNORDERED_SET_HM_H

#include <ext/pb_ds/assoc_container.hpp>
#include <tr1/unordered_set>

namespace piranha
{
  template <class Element, class Hasher, class Eq, class Allocator, bool StoreHash>
    class mult_hash
  {
      typedef std::tr1::unordered_set<Element,Hasher,Eq,Allocator,StoreHash> container_type;
    public:
      typedef typename container_type::const_iterator iterator;
      typedef typename container_type::const_iterator point_iterator;
      mult_hash()
      {
          private_container_.max_load_factor(.5);
      }
      mult_hash(const size_t &s):private_container_(s)
      {
          private_container_.max_load_factor(.5);
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
    class c_mult_hash_
  {
      typedef pb_ds::cc_hash_table<
        Element,
        pb_ds::null_mapped_type,
        Hasher,
        Eq,
        pb_ds::direct_mask_range_hashing<size_t>,
        pb_ds::hash_standard_resize_policy<
          pb_ds::hash_exponential_size_policy<size_t>,
          pb_ds::hash_load_check_resize_trigger<false,size_t>,
          true,
          size_t
        >,
        true,
        Allocator
      > container_type;
    public:
      typedef typename container_type::const_iterator iterator;
      typedef typename container_type::point_iterator point_iterator;
      c_mult_hash_()
        {}
      c_mult_hash_(const size_t &)
        {}
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
    typedef c_mult_hash_<Element,Hasher,Eq,Allocator,StoreHash> type;
  };
}

#endif
