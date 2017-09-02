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

#ifndef PIRANHA_CODED_HASH_TABLE_H
#define PIRANHA_CODED_HASH_TABLE_H

#include <algorithm>
#include <boost/cstdint.hpp>
#include <boost/functional/hash.hpp>
#include <boost/integer_traits.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <exception>
#include <utility> // For std::pair.
#include <vector>

#include "config.h"
#include "exceptions.h"
#include "settings.h"

namespace piranha
{

// TODO: make sure the size of the vector never goes past MaxFastInt, since we are converting size_type -> MaxFastInt in the hash coded
// multiplier, when determining the memory position of the key.

/// Hash table for highly-sparse coded series multiplication.
template <class Cf, class Code, class Allocator>
class coded_hash_table
{
        public:
                /// The Code is the key type for the hash table.
                typedef Code KeyType;
                /// Value type.
                typedef std::pair<Cf,Code> value_type;
                /// Bucket.
                /**
                 * Bucket size N is fixed at compile time and must be in the ]0,256] range.
                 */
                template <int N>
                class bucket
                {
                        public:
                                /// Bucket size type.
                                typedef boost::uint8_t size_type;
                static_assert(N > 0 && N <= boost::integer_traits<size_type>::const_max, "Invalid bucket size.");
                                /// Default constructor.
                                /**
                                 * Will initialise the code of each coefficient to -1.
                                 */
                                bucket() 
                                {
                                        for (size_type i = 0; i < N; ++i) {
                                                t[i].second = -1;
                                        }
                                }
                                /// Array of coefficient-code pairs.
                                value_type t[N];
                };
                /// Size policy.
                enum size_policy {
                        /// Power of two sizes.
                        pow2    = 0,
                        /// Prime sizes.
                        prime   = 1
                };
                // Configuration options.
                /// Bucket size.
                static const std::size_t bucket_size            = 5;
                /// Minimum size index.
                static const std::size_t min_size_index         = 0;
                /// Maximum number of buckets probed outside the native bucket.
                static const std::size_t probe_size             = 100;
                // Configuration options end here.
                static const std::size_t sizes_size =
#ifdef _PIRANHA_64BIT
                        40;
#else
                        32;
#endif
                /// List of available sizes for the hash table.
                /**
                 * The number of available sizes depends on the architecture.
                 */
                static const std::size_t sizes[2][sizes_size];
                /// Alias for the bucket type.
                typedef bucket<bucket_size> bucket_type;
                /// Bucket size type.
                typedef typename bucket_type::size_type bucket_size_type;
                /// Internal container type.
                typedef std::vector<bucket_type,typename Allocator::template rebind<bucket_type>::other> container_type;
                /// Size type.
                typedef typename container_type::size_type size_type;
                /// Iterator class.
                class iterator: public boost::iterator_facade<iterator,value_type,boost::forward_traversal_tag>
                {
                                friend class boost::iterator_core_access;
                                friend class coded_hash_table;
                                /// Constructor from coded_hash_table.
                                /**
                                 * The iterator will point to the first element or to the end of the table, if
                                 * the table is empty.
                                 *
                                 * @param[in] p pointer to the hash table.
                                 */
                                iterator(coded_hash_table *p): m_ht(p), m_vector_index(0), m_bucket_index(0)
                                {
                                        // If the first slot is not taken, find the next one.
                                        if (m_ht->m_container[m_vector_index].t[m_bucket_index].second < 0) {
                                                increment();
                                        }
                                }
                                /// Constructor from coded_hash_table, and vector-bucket indices.
                                /**
                                 * No checks are performed on the values supplied to the constructor.
                                 *
                                 * @param[in] p pointer to the hash table.
                                 * @param[in] vi vector index.
                                 * @param[in] bi bucket index.
                                 */
                                iterator(coded_hash_table *p, const size_type &vi, const bucket_size_type &bi):
                                        m_ht(p), m_vector_index(vi), m_bucket_index(bi) {}
                        public:
                                iterator():m_ht(0), m_vector_index(0), m_bucket_index(0) {}
                        private:
                                value_type &dereference() const
                                {
                                        PIRANHA_ASSERT(m_ht);
                                        PIRANHA_ASSERT(m_vector_index < m_ht->m_container.size());
                                        PIRANHA_ASSERT(m_bucket_index < bucket_size);
                                        PIRANHA_ASSERT(m_ht->m_container[m_vector_index].t[m_bucket_index].second >= 0);
                                        return m_ht->m_container[m_vector_index].t[m_bucket_index];
                                }
                                bool equal(const iterator &it2) const
                                {
                                        return (m_ht == it2.m_ht && m_vector_index == it2.m_vector_index &&
                                                m_bucket_index == it2.m_bucket_index);
                                }
                                void increment()
                                {
                                        const size_type vector_size = m_ht->m_container.size();
                                        while (true) {
                                                // Go to the next bucket if we are at the last element of the current one.
                                                if (m_bucket_index == bucket_size - 1) {
                                                        m_bucket_index = 0;
                                                        ++m_vector_index;
                                                } else {
                                                        // Otherwise just go to the next element of the bucket.
                                                        ++m_bucket_index;
                                                }
                                                // If we went past the vector size or if we found the next element, break out and return.
                                                if (m_vector_index == vector_size || m_ht->m_container[m_vector_index].t[m_bucket_index].second >= 0) {
                                                        break;
                                                }
                                        }
                                }
                        private:
                                coded_hash_table        *m_ht;
                                size_type               m_vector_index;
                                bucket_size_type        m_bucket_index;
                };
                /// Default constructor.
                /**
                 * Sets the size policy to pow2.
                 */
                coded_hash_table(): m_size_policy(pow2),m_size_index(min_size_index),m_length(0),
                        m_container(boost::numeric_cast<size_type>(sizes[m_size_policy][m_size_index]))
                {}
                // NOTE: remove?
                /// Constructor with size hint.
                /**
                 * Sets the size policy to pow2. The hash table size is set to a value close to the provided size hint.
                 *
                 * @param[in] size_hint size hint for the hash table.
                 */
                coded_hash_table(const std::size_t &size_hint):
                        m_size_policy(pow2),m_size_index(find_upper_size_index(size_hint / bucket_size + 1)),m_length(0),
                        m_container(boost::numeric_cast<size_type>(sizes[m_size_policy][m_size_index]))
                {}
                /// Destructor.
                ~coded_hash_table()
                {
                        PIRANHA_DEBUG(std::cout << "On destruction, the vector size of coded_hash_table was: "
                                << sizes[m_size_policy][m_size_index] << '\n');
                        PIRANHA_DEBUG(std::cout << "On destruction, the load factor was: "
                                << double(m_length) / (sizes[m_size_policy][m_size_index] * bucket_size) << '\n');
// GUT: commented . throws inside a destructor
//                      PIRANHA_ASSERT(sizes[m_size_policy][m_size_index] == m_container.size());
                }
                /// Return an iterator to the first element of the table.
                /**
                 * If the table is empty, the end() iterator is returned.
                 *
                 * @return iterator to the beginning of the table.
                 */
                iterator begin()
                {
                        return iterator(this);
                }
                /// Return an iterator to the end of the table.
                /**
                 *
                 * @return iterator to the end of the series.
                 */
                iterator end()
                {
                        return iterator(this,m_container.size(),0);
                }
                /// Locate element based on key.
                /**
                 * If the element is found, then return (true,position). Otherwise return (false, first empty slot). Note that in the latter
                 * case, the first empty slot can be the table's end.
                 *
                 * @param[in] key value key.
                 *
                 * @return pair containing the result of the operation.
                 */
                std::pair<bool,iterator> find(const KeyType &key)
                {
                        const size_type vector_size = m_container.size(),
                                vector_pos = get_position(boost::hash<KeyType>()(key),vector_size,m_size_policy);
                        PIRANHA_ASSERT(vector_pos < vector_size);
                        const bucket_type &bucket = m_container[vector_pos];
                        // Now examine all elements in the bucket.
                        for (bucket_size_type i = 0; i < bucket_size; ++i) {
                                // If the slot in the bucket is not taken (which means there are no more elements to examine),
                                // it means that key was not found.
                                if (bucket.t[i].second < 0) {
                                        return std::make_pair(false,iterator(this, vector_pos, i));
                                }
                                if (bucket.t[i].second == key) {
                                        // If we found an occupied bucket slot, examine the key to see whether it matches or not
                                        // with t's.
                                        return std::make_pair(true,iterator(this, vector_pos, i));
                                }
                                // No bucket end and no match, continue to the next bucket element.
                        }
                        // We examined all the elements in the (full) destination bucket, found no match. Start linear probing.
                        // Start looking from the next bucket.
                        size_type new_vector_pos = vector_pos + 1;
                        for (size_type n = 0; n < probe_size; ++n, ++new_vector_pos) {
                                // Break out if we are at the end of the table.
                                if (new_vector_pos == vector_size) {
                                        break;
                                }
                                // Look into the bucket, as above.
                                const bucket_type &bucket = m_container[new_vector_pos];
                                for (bucket_size_type i = 0; i < bucket_size; ++i) {
                                        if (bucket.t[i].second < 0) {
                                                return std::make_pair(false,iterator(this, new_vector_pos, i));
                                        }
                                        if (bucket.t[i].second == key) {
                                                return std::make_pair(true,iterator(this, new_vector_pos, i));
                                        }
                                }
                        }
                        // Either we exhausted the linear probe sequence or we are at the end of the table. Return (false,end()).
                        return std::make_pair(false,end());
                }
                /// Insert new value into the hash table at the position specified by iterator.
                /**
                 * Iterator must point either to an empty slot or to the end of the table. Value is assumed not to be already present in the table.
                 *
                 * @param[in] v value_type to be inserted.
                 * @param[in] it insertion point.
                 */
                void insert_new(const value_type &v, const iterator &it)
                {
                        if (!attempt_insertion(v,it)) {
                                iterator tmp(it);
                                do {
                                        PIRANHA_DEBUG(std::cout << "Started resizing coded series hash table." << '\n');
                                        increase_size();
                                        tmp = find(v.second).second;
                                        PIRANHA_DEBUG(std::cout << "Resized coded series hash table." << '\n');
                                } while (!attempt_insertion(v,tmp));
                        }
                }
                /// Return the number of elements stored in the hash table.
                /**
                 * @return the number of items in the hash table.
                 */
                size_type size() const
                {
                        return m_length;
                }
                /// Get memory position in which key would be inserted, relative to the table's starting point.
                size_type get_memory_position(const KeyType &key) const
                {
                        return get_position(boost::hash<KeyType>()(key),m_container.size(),m_size_policy);
                }
                /// Return the size of the vector representing the hash table internally.
                size_type get_vector_size() const
                {
                        return m_container.size();
                }
        private:
                // Convert hash value into position.
                static size_type get_position(const std::size_t &hash, const size_type &size, const size_policy &sp)
                {
                        switch (sp) {
                                case pow2:
                                        return (hash & (size - 1));
                                case prime:
                                        return (hash % size);
                                default:
                                        PIRANHA_ASSERT(false);
                                        return 0;
                        }
                }
                // TODO: remove?
                // Find hash table size from hint.
                size_type find_upper_size_index(const std::size_t &size) const
                {
                        for (std::size_t retval = 0; retval < sizes_size; ++retval) {
                                if (sizes[m_size_policy][retval] >= size) {
                                        return std::max<std::size_t>(min_size_index,retval);
                                }
                        }
                        return min_size_index;
                }
                // Attempt insertion and return the outcome. Iterator must point either to an empty slot or to the end of the table.
                bool attempt_insertion(const value_type &v, const iterator &it)
                {
                        const size_type vector_size = m_container.size(), vector_index = it.m_vector_index;
                        // If the iterator is the end, we cannot insert.
                        if (vector_index == vector_size)
                        {
                                return false;
                        }
                        const bucket_size_type bucket_index = it.m_bucket_index;
                        PIRANHA_ASSERT(bucket_index < bucket_size);
                        PIRANHA_ASSERT(vector_index < vector_size);
                        bucket_type &bucket = m_container[vector_index];
                        // Make sure the slot is not already taken.
                        PIRANHA_ASSERT(bucket.t[bucket_index].second < 0);
                        bucket.t[bucket_index] = v;
                        ++m_length;
                        return true;
                }
                // Insertion routine that won't check for equal key. Used during resizing. Outcome is reported.
                bool unchecked_insertion(const value_type &v)
                {
                        const size_type vector_size = m_container.size(),
                                vector_pos = get_position(boost::hash<KeyType>()(v.second),vector_size,m_size_policy);
                        PIRANHA_ASSERT(vector_pos < vector_size);
                        bucket_type &bucket = m_container[vector_pos];
                        // Now check for an available slot in the bucket.
                        for (bucket_size_type i = 0; i < bucket_size; ++i) {
                                // If the slot in the bucket is not taken we can place the key here.
                                if (bucket.t[i].second < 0) {
                                        bucket.t[i] = v;
                                        ++m_length;
                                        return true;
                                }
                        }
                        // Try insert into the next bucket(s).
                        size_type new_vector_pos = vector_pos + 1;
                        for (size_type n = 0; n < probe_size; ++n, ++new_vector_pos) {
                                if (new_vector_pos == vector_size) {
                                        break;
                                }
                                bucket_type &bucket = m_container[new_vector_pos];
                                for (bucket_size_type i = 0; i < bucket_size; ++i) {
                                        // If slot is not taken, use it.
                                        if (bucket.t[i].second < 0) {
                                                bucket.t[i] = v;
                                                ++m_length;
                                                return true;
                                        }
                                }
                        }
                        // We were not able to insert.
                        return false;
                }
                // Increase size of the container to the next size.
                void increase_size()
                {
                        const double load_factor = static_cast<double>(m_length) / (static_cast<double>(m_container.size()) * bucket_size);
                        PIRANHA_DEBUG(std::cout << "Increase size requested at load factor: " << load_factor << '\n');
                        // If load factor is too small and we are on pow2 sizes,
                        // we want to switch to prime sizes.
                        static const double min_load_factor = .1;
                        size_policy new_size_policy = m_size_policy;
                        if (m_size_policy == pow2 && load_factor < min_load_factor) {
                                PIRANHA_DEBUG(std::cout << "Load factor too low in pow2 sizes, switching to prime sizes.\n";)
                                new_size_policy = prime;
                        }
                        coded_hash_table new_ht;
                        new_ht.m_size_policy = new_size_policy;
                        new_ht.m_size_index = m_size_index + 2;
                        new_ht.m_length = 0;
                        new_ht.m_container.resize(boost::numeric_cast<size_type>(sizes[new_ht.m_size_policy][new_ht.m_size_index]));
                        // Cache quantities.
                        const iterator it_i = begin(), it_f = end();
                        iterator it = it_i;
                        size_type count = 0;
                        while (it != it_f) {
                                if (!new_ht.unchecked_insertion(*it)) {
                                        // NOTICE: here maybe we can use swapping instead of copying. The only problem is that
                                        // resizing can fail. In that case, we should swap back everything, if possible, and re-attempt
                                        // the resize with a bigger value.
                                        PIRANHA_DEBUG(std::cout << "Hash table resize triggered during resize." << '\n');
                                        // If we are able to rebuild less than a fraction of the initial hash table
                                        // and we are working in pow2 sizes, switch to prime sizes.
                                        static const double rebuild_thresh = .5;
                                        if (new_ht.m_size_policy == pow2 && static_cast<double>(count) / m_length < rebuild_thresh) {
                                                PIRANHA_DEBUG(std::cout << "Rebuilding failed too early, switching to prime sizes.\n");
                                                new_ht.m_size_policy = prime;
                                        }
                                        ++new_ht.m_size_index;
                                        if (new_ht.m_size_index == sizes_size) {
                                                // We ran out of possible sizes.
                                                PIRANHA_THROW(std::overflow_error,"hash table size overflow");
                                        }
                                        new_ht.m_length = 0;
                                        new_ht.m_container.clear();
                                        new_ht.m_container.resize(
                                                boost::numeric_cast<size_type>(sizes[new_ht.m_size_policy][new_ht.m_size_index]));
                                        it = it_i;
                                        count = 0;
                                } else {
                                        ++it;
                                        ++count;
                                }
                        }
                        // Rebuilding complete, swap data members.
                        m_container.swap(new_ht.m_container);
                        std::swap(m_size_index,new_ht.m_size_index);
                        std::swap(m_length,new_ht.m_length);
                        std::swap(m_size_policy,new_ht.m_size_policy);
                }
        private:
                size_policy             m_size_policy;
                std::size_t             m_size_index;
                size_type               m_length;
                container_type          m_container;
};

template <class Cf, class Code, class Allocator>
const std::size_t coded_hash_table<Cf,Code,Allocator>::min_size_index;

template <class Cf, class Code, class Allocator>
const std::size_t coded_hash_table<Cf,Code,Allocator>::bucket_size;

template <class Cf, class Code, class Allocator>
const std::size_t coded_hash_table<Cf,Code,Allocator>::sizes[2][coded_hash_table<Cf,Code,Allocator>::sizes_size] = { {
        1ull,
        2ull,
        4ull,
        8ull,
        16ull,
        32ull,
        64ull,
        128ull,
        256ull,
        512ull,
        1024ull,
        2048ull,
        4096ull,
        8192ull,
        16384ull,
        32768ull,
        65536ull,
        131072ull,
        262144ull,
        524288ull,
        1048576ull,
        2097152ull,
        4194304ull,
        8388608ull,
        16777216ull,
        33554432ull,
        67108864ull,
        134217728ull,
        268435456ull,
        536870912ull,
        1073741824ull,
        2147483648ull
#ifdef _PIRANHA_64BIT
        ,
        4294967296ull,
        8589934592ull,
        17179869184ull,
        34359738368ull,
        68719476736ull,
        137438953472ull,
        274877906944ull,
        549755813888ull
#endif
}, {
        1ull,
        3ull,
        5ull,
        11ull,
        23ull,
        53ull,
        97ull,
        193ull,
        389ull,
        769ull,
        1543ull,
        3079ull,
        6151ull,
        12289ull,
        24593ull,
        49157ull,
        98317ull,
        196613ull,
        393241ull,
        786433ull,
        1572869ull,
        3145739ull,
        6291469ull,
        12582917ull,
        25165843ull,
        50331653ull,
        100663319ull,
        201326611ull,
        402653189ull,
        805306457ull,
        1610612741ull,
        3221225473ull
#ifdef _PIRANHA_64BIT
        ,
        6442450939ull,
        12884901893ull,
        25769803799ull,
        51539607551ull,
        103079215111ull,
        206158430209ull,
        412316860441ull,
        824633720831ull
#endif
} };
}

#endif
