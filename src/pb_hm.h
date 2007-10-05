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

#ifndef PIRANHA_PB_HM_H
#define PIRANHA_PB_HM_H

#if GCC_VERSION < 402000
#include <ext/pb_assoc/assoc_cntnr.hpp>
#else
#include <ext/pb_ds/assoc_container.hpp>
#endif

namespace piranha
{
  template <class Element, class Hasher, class Eq, class Allocator, bool StoreHash>
    class mult_hash
  {
#if GCC_VERSION < 402000
      typedef pb_assoc::cc_hash_assoc_container<
        Element,
        pb_assoc::null_mapped_type,
        Hasher,
        std::equal_to<Element>,
        pb_assoc::direct_mask_range_hashing<size_t>,
        pb_assoc::hash_standard_resize_policy<
          pb_assoc::hash_exponential_size_policy<size_t>,
          pb_assoc::hash_load_check_resize_trigger<false,size_t>,
          true,
          size_t
        >,
        true,
        Allocator
      > container_type;
#else
      typedef pb_ds::cc_hash_table<
        Element,
        pb_ds::null_mapped_type,
        Hasher,
        std::equal_to<Element>,
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
#endif
    public:
      typedef typename container_type::const_iterator iterator;
      typedef typename container_type::point_iterator point_iterator;
      mult_hash()
        {}
      mult_hash(const size_t &)
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
}

#endif
