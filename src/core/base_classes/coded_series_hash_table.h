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

#ifndef PIRANHA_CODED_SERIES_HASH_TABLE_H
#define PIRANHA_CODED_SERIES_HASH_TABLE_H

#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>
#include <vector>

#include "../p_assert.h"
#include "../settings.h"

namespace piranha
{
  // Forward declaration.
  template <class, class, int>
    class coded_series_hash_table;

  template <class Cf, class Ckey, int N>
    class coded_term_bucket
  {
      BOOST_STATIC_ASSERT(N > 0);
      friend class coded_series_hash_table<Cf,Ckey,N>;
    public:
      static const size_t size = (size_t)N;
      struct term
      {
        mutable Cf  m_cf;
        Ckey        m_ckey;
      };
      coded_term_bucket() {init();}
    private:
      void init()
      {
        for (size_t i = 0; i < size; ++i)
        {
          m_flags[i] = false;
        }
      }
    private:
      term  m_terms[N];
      bool  m_flags[N];
  };

  template <class Cf, class Ckey, int N>
    class coded_series_hash_table
  {
      typedef coded_term_bucket<Cf,Ckey,N> bucket_type;
      typedef std::vector<bucket_type> container_type;
      static const size_t initial_vector_size = 32;
      BOOST_STATIC_ASSERT(initial_vector_size > 0);
      BOOST_STATIC_ASSERT((initial_vector_size & (initial_vector_size - 1)) == 0);
      static const size_t bucket_size = bucket_type::size;
    public:
      typedef typename bucket_type::term term_type;
      class iterator
      {
        public:
          iterator(const coded_series_hash_table *p):m_ht(p),m_vector_index(0),m_bucket_index(0)
          {
            // If the first slot is not taken, find the next one.
            if (!m_ht->m_container[m_vector_index].m_flags[m_bucket_index])
            {
              next();
            }
          }
          iterator(const coded_series_hash_table *p,const size_t &vi, const size_t &bi):
            m_ht(p),m_vector_index(vi),m_bucket_index(bi) {}
          iterator &operator++()
          {
            next();
            return *this;
          }
          const term_type &operator*() const
          {
            p_assert(m_vector_index < m_ht->m_container.size());
            p_assert(m_bucket_index < bucket_size);
            p_assert(m_ht->m_container[m_vector_index].m_flags[m_bucket_index]);
            return m_ht->m_container[m_vector_index].m_terms[m_bucket_index];
          }
          const term_type *operator->() const
          {
            p_assert(m_vector_index < m_ht->m_container.size());
            p_assert(m_bucket_index < bucket_size);
            p_assert(m_ht->m_container[m_vector_index].m_flags[m_bucket_index]);
            return &m_ht->m_container[m_vector_index].m_terms[m_bucket_index];
          }
          bool operator==(const iterator &it2) const
          {
            return (m_ht == it2.m_ht and m_vector_index == it2.m_vector_index and
              m_bucket_index == it2.m_bucket_index);
          }
          bool operator!=(const iterator &it2) const
          {
            return (m_ht != it2.m_ht or m_vector_index != it2.m_vector_index or
              m_bucket_index != it2.m_bucket_index);
          }
        private:
          void next()
          {
            const size_t vector_size = m_ht->m_container.size();
            while (true)
            {
              // Go to the next bucket if we are at the last element of the current one.
              if (m_bucket_index == bucket_size - 1)
              {
                m_bucket_index = 0;
                ++m_vector_index;
              }
              // Otherwise just go to the next element of the bucket.
              else
              {
                ++m_bucket_index;
              }
              // If we went past the vector size or if we found the next element, break out and return.
              if (m_vector_index == vector_size or m_ht->m_container[m_vector_index].m_flags[m_bucket_index])
              {
                break;
              }
            }
          }
        private:
          const coded_series_hash_table *m_ht;
          size_t                        m_vector_index;
          size_t                        m_bucket_index;
      };
      coded_series_hash_table():m_container(initial_vector_size) {}
      iterator begin() const
      {
        return iterator(this);
      }
      iterator end() const
      {
        return iterator(this,m_container.size(),0);
      }
      iterator find(const term_type &t) const
      {
        // TODO: maybe it is more efficient to call hash_value directly here?
        boost::hash<Ckey> ckey_hash;
        const size_t h = ckey_hash(t.m_ckey), vector_size = m_container.size();
        // Check that the size is a power of two.
        p_assert(vector_size > 0 and (vector_size & (vector_size - 1)) == 0);
        // Find h % vector_size using logical AND, since we are working with powers of two.
        const size_t vector_pos = h & (vector_size - 1);
        p_assert(vector_pos < vector_size);
        // Now examine all elements in the bucket.
        for (size_t i = 0; i < bucket_size; ++i)
        {
          // If the slot in the bucket is not taken (which means there are no more elements to examine),
          // it means that t was not found.
          if (!m_container[vector_pos].m_flags[i])
          {
            return end();
          }
          // If we found an occupied bucket slot, examine the key to see whether it matches or not with t's.
          // If it does not match, let's move to the next bucket element.
          else if (m_container[vector_pos].m_terms[i].m_ckey == t.m_ckey)
          {
            return iterator(this,vector_pos,i);
          }
        }
        // All the elements of the bucket were taken, we examined them but found no match.
        return end();
      }
      void insert(const term_type &t)
      {
        while (!attempt_insertion(t))
        {
          // Increase size until insertion succeeds.
          increase_size();
          __PDEBUG(std::cout << "Resized coded series hash table." << '\n');
        }
      }
    private:
      bool attempt_insertion(const term_type &t)
      {
        // TODO: for the future, when we are sure this structure works well: would it be more
        // efficient if insert() automatically took charge of updating the value if a term with equal key
        // is found? Just like base_series does... Maybe not, there seems to be no difference.
        // TODO: maybe it is more efficient to call hash_value directly here?
        boost::hash<Ckey> ckey_hash;
        const size_t h = ckey_hash(t.m_ckey), vector_size = m_container.size();
        // Check that the size is a power of two.
        p_assert(vector_size > 0 and (vector_size & (vector_size - 1)) == 0);
        // Find h % vector_size using logical AND, since we are working with powers of two.
        const size_t vector_pos = (h & (vector_size - 1));
        p_assert(vector_pos < vector_size);
        // Now examine all elements in the bucket.
        for (size_t i = 0; i < bucket_size; ++i)
        {
          // If the slot in the bucket is not taken (which means there are no more elements to examine),
          // it means that we found a place for t.
          if (!m_container[vector_pos].m_flags[i])
          {
            // Set the flag to occupied.
            m_container[vector_pos].m_flags[i] = true;
            // Let's copy t over the correct position.
            // TODO: maybe here a swap() can be faster.
            m_container[vector_pos].m_terms[i] = t;
            return true;
          }
        }
        return false;
      }
      // Increase size of the container by a factor of a natural power of 2.
      void increase_size()
      {
        size_t power = 1;
        coded_series_hash_table new_ht;
        new_ht.m_container.resize(m_container.size() << power);
        const iterator it_i = begin(), it_f = end();
        for (iterator it = it_i; it != it_f; ++it)
        {
          if (!new_ht.attempt_insertion(*it))
          {
            ++power;
            new_ht.m_container.clear();
            new_ht.m_container.resize(new_ht.m_container.size() << power);
            it = it_i;
          __PDEBUG(std::cout << "Hash table resize triggered during resize." << '\n');
          }
        }
        new_ht.m_container.swap(m_container);
      }
    public:
      container_type  m_container;
  };
}

#endif
