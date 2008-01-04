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

#ifndef PIRANHA_SYMBOL_LIMITER_H
#define PIRANHA_SYMBOL_LIMITER_H

#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <deque>

#include "psymbol.h"

namespace piranha
{
  class symbol_limiter
  {
    public:
      static void set_limit(psym_p, uint16);
      static void set_limit(const std::string &, uint16);
      static void clear();
      static void put();
    private:
      symbol_limiter()
        {}
// Boilerplate for the multiindex container.
      struct limit_element
      {
        limit_element(psym_p s, uint16 n):symbol(s),limit(n)
          {}
        psym_p      symbol;
        uint16   limit;
      };
// We can use this structure as comparison functor, whereas below we hash according to name,
// because in managing symbols we make sure that there are no symbols with same name. Hence
// symbol pointer and name are effectively the same thing.
      struct eq_psym_p
      {
        bool operator()(psym_p p1, psym_p p2) const
        {
          return (p1==p2);
        }
      };
      struct hash_psym_p
      {
        size_t operator()(psym_p p) const
        {
          return boost::hash<std::string>()(p->name());
        }
      };
      struct limit_modifier
      {
        limit_modifier(uint16 n):new_n_(n)
          {}
        ~limit_modifier()
          {}
        void operator()(limit_element &l)
        {
          l.limit = new_n_;
        }
        uint16   new_n_;
      };
      typedef boost::multi_index_container < limit_element,
        boost::multi_index::indexed_by <
        boost::multi_index::hashed_unique <
        boost::multi_index::member<limit_element,psym_p,&limit_element::symbol>,hash_psym_p,eq_psym_p
        >
        >
        >
        limits_map;
      typedef limits_map::nth_index<0>::type::iterator map_iterator;
    public:
/// Limits' indices.
/**
 * When constructed from a vector of pointers to symbols, it will build a vector specifying limits
 * in terms of the symbols' positions in the input vector.
 */
      class index_limit
      {
        typedef std::deque<boost::tuple<size_t,uint16> > vec_expo_index_limit;
        public:
          index_limit(const vector_psym_p &v):private_min_limit_(0),private_limits_()
          {
            if (lmap_.empty())
            {
              return;
            }
// Initial limit is the first element of limit map.
            private_min_limit_=lmap_.begin()->limit;
            const size_t w=v.size();
            map_iterator mit;
            const map_iterator mit_f=lmap_.end();
            boost::tuple<size_t,uint16> insert_tuple;
            for (size_t i=0;i<w;++i)
            {
              mit=find_expo_limit(v[i]);
              if (mit!=mit_f)
              {
                insert_tuple.get<0>()=i;
                insert_tuple.get<1>()=mit->limit;
                private_limits_.push_back(insert_tuple);
                if (insert_tuple.get<1>() < private_min_limit_)
                {
                  private_min_limit_ = insert_tuple.get<1>();
                }
              }
            }
          }
          const boost::tuple<size_t,uint16> &operator[](const size_t &n) const
          {
            return private_limits_[n];
          }
          size_t size() const
          {
            return private_limits_.size();
          }
          const uint16 &g_min_expo() const
          {
            return private_min_limit_;
          }
        private:
// Make default ctor private, just to make sure it is not called.
          index_limit() {}
// Data members.
          uint16             private_min_limit_;
          vec_expo_index_limit  private_limits_;
      };
    private:
/// Check whether a limit has already been set for a symbol.
      static map_iterator find_expo_limit(psym_p it)
      {
        return lmap_.find(it);
      }
// Data members.
    private:
      static limits_map   lmap_;
  };

// Raise index_limit to piranha namespace.
  typedef symbol_limiter::index_limit index_limit;
}
#endif
