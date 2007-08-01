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

#ifndef PIRANHA_SYMBOL_LIMITER_H
#define PIRANHA_SYMBOL_LIMITER_H

#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>

//TODO: use multiindex container, so that we can drop the dependency on TR1.

#include "psymbol.h"

namespace piranha
{
  typedef std::vector<boost::tuple<size_t,int> > vec_expo_index_limit;

  class symbol_limiter
  {
    public:
      static void set_limit(psym_p, int);
      static void set_limit(const std::string &, int);
      static void clear();
      static void get_limits_index(const vector_psym_p &v, vec_expo_index_limit &retval)
      {
        p_assert(retval.size()==0);
        const size_t w=v.size();
// Reserve space so that we potentially save resizes.
// TODO: check the real impact of this by commenting and benchmarking.
        retval.reserve(w);
        map_iterator mit;
        const map_iterator mit_f=lmap_.end();
        boost::tuple<size_t,int> insert_tuple;
        for (size_t i=0;i<w;++i)
          {
            mit=find_expo_limit(v[i]);
            if (mit!=mit_f)
              {
                insert_tuple.get<0>()=i;
                insert_tuple.get<1>()=mit->limit;
                retval.push_back(insert_tuple);
              }
          }
      }
      static void put();
    private:
      symbol_limiter()
        {}
// Boilerplate for the multiindex container.
      struct limit_element
      {
        limit_element(psym_p s, int n):symbol(s),limit(n)
        {}
        psym_p  symbol;
        int     limit;
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
        limit_modifier(int n):new_n_(n)
          {}
        ~limit_modifier()
          {}
        void operator()(limit_element &l)
        {
          l.limit = new_n_;
        }
        int new_n_;
      };

    public:
      typedef boost::multi_index_container < limit_element,
        boost::multi_index::indexed_by <
          boost::multi_index::hashed_unique <
            boost::multi_index::member<limit_element,psym_p,&limit_element::symbol>,hash_psym_p,eq_psym_p
          >
        >
      >
      limits_map;
      //typedef limits_map::iterator map_iterator;
      typedef limits_map::nth_index<0>::type::iterator map_iterator;
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
}

#endif
