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

#ifndef PIRANHA_CHAINED_HASH_SET_H
#define PIRANHA_CHAINED_HASH_SET_H

#include <algorithm> // For std::swap.
#include <boost/functional/hash.hpp> // Default hasher.
#include <functional> // std::equal_to.
#include <memory> // std::allocator.
#include <vector>

#include "common_functors.h"
#include "config.h"
#include "memory.h"
#include "p_assert.h"
#include "settings.h"

namespace piranha
{
	template <class T, class HashFunction = boost::hash<T>, class EqualKey = std::equal_to<T>,
	class Allocator = std::allocator<char>, class SwapFunction = std_swap<T> >
	class chained_hash_set
	{
			typedef std::vector<T,counting_allocator<T,Allocator> > bucket_type;
			typedef std::vector<bucket_type,counting_allocator<T,Allocator> > vector_type;
			// Configuration options.
			static const size_t min_size_index =	1;
			static const size_t shrink_trigger = 	10;
			// Configuration options stop here.
			static const size_t sizes_size =
#ifdef _PIRANHA_64BIT
				40;
#else
				32;
#endif
			static const size_t sizes[sizes_size];
			p_static_check(shrink_trigger > 1, "Shrink trigger must be strictly greater than 1.");
			p_static_check(min_size_index > 0, "Min size index must be at least 1.");
		public:
			// Public typedefs.
			typedef T value_type;
			typedef T key_type;
			typedef HashFunction hasher;
			typedef EqualKey key_equal;
			typedef SwapFunction swapper;
			typedef T * pointer;
			typedef const T * const_pointer;
			typedef T & reference;
			typedef const T & const_reference;
			typedef size_t size_type;
			typedef ptrdiff_t difference_type;
			class iterator
			{
					friend class chained_hash_set;
					explicit iterator(const chained_hash_set *ht): m_ht(ht), m_vector_index(0), m_bucket_index(0) {
						// If the first bucket is empty taken, find the first element of the hash set.
						if (m_ht->m_container[m_vector_index].size() == 0) {
							next();
						}
					}
					explicit iterator(const chained_hash_set *ht, const size_t &vi, const size_t &bi):
							m_ht(ht), m_vector_index(vi), m_bucket_index(bi) {}
				public:
					iterator &operator++() {
						next();
						return *this;
					}
					const key_type &operator*() const {
						p_assert(m_vector_index < m_ht->m_container.size());
						p_assert(m_bucket_index < m_ht->m_container[m_vector_index].size());
						return m_ht->m_container[m_vector_index][m_bucket_index];
					}
					const key_type *operator->() const {
						p_assert(m_vector_index < m_ht->m_container.size());
						p_assert(m_bucket_index < m_ht->m_container[m_vector_index].size());
						return &m_ht->m_container[m_vector_index][m_bucket_index];
					}
					bool operator==(const iterator &it2) const {
						return (m_ht == it2.m_ht && m_vector_index == it2.m_vector_index &&
								m_bucket_index == it2.m_bucket_index);
					}
					bool operator!=(const iterator &it2) const {
						return !(*this == it2);
					}
				private:
					void next() {
						const size_t vector_size = m_ht->m_container.size();
						++m_bucket_index;
						// If this was the last element of the bucket or we incremented the index in a
						// empty bucket, we need to move to the next bucket (which we do not know if it
						// is empty or not, yet).
						p_assert(m_vector_index < vector_size);
						if (m_bucket_index >= m_ht->m_container[m_vector_index].size()) {
							m_bucket_index = 0;
							++m_vector_index;
						}
						// Let's find the first non-empty bucket, if needed, and stop at the end of the hash set.
						while (m_vector_index < vector_size && m_ht->m_container[m_vector_index].size() == 0) {
							m_bucket_index = 0;
							++m_vector_index;
						}
					}
				private:
					chained_hash_set const	*m_ht;
					size_t					m_vector_index;
					size_t					m_bucket_index;
			};
			typedef iterator const_iterator;
			chained_hash_set(): m_sizes_index(min_size_index), m_length(0), m_container(sizes[m_sizes_index])
			{}
			chained_hash_set(const size_t &hint):
				m_sizes_index(find_upper_size_index(hint/settings::load_factor())),
				m_length(0),
				m_container(sizes[m_sizes_index])
			{}
			~chained_hash_set() {
				__PDEBUG(std::cout << "On destruction, the vector size of chained_hash_set was: "
								   << m_container.size() << '\n');
			}
			iterator begin() const {
				return iterator(this);
			}
			iterator end() const {
				return iterator(this, sizes[m_sizes_index], 0);
			}
			iterator find(const key_type &key) const {
				const size_t vector_pos = get_position(hasher()(key),sizes[m_sizes_index]);
				p_assert(vector_pos < sizes[m_sizes_index]);
				// Now examine all elements in the bucket.
				const size_t bucket_size = m_container[vector_pos].size();
				for (size_t i = 0; i < bucket_size; ++i) {
					if (m_container[vector_pos][i] == key) {
						return iterator(this, vector_pos, i);
					}
				}
				// No matching element was found in the bucket.
				return end();
			}
			iterator insert(const key_type &key) {
				p_assert(settings::load_factor() > 0);
				if (static_cast<double>(m_length + 1) > settings::load_factor() * sizes[m_sizes_index]) {
					grow();
				}
				const size_t vector_pos = get_position(hasher()(key),sizes[m_sizes_index]);
				p_assert(vector_pos < sizes[m_sizes_index]);
				// Append the key at the end of the bucket.
				// TODO: think about placing it in front.
				// OPTIONS: use std::deque or swap with front.
				m_container[vector_pos].push_back(key);
				++m_length;
				return iterator(this,vector_pos,m_container[vector_pos].size() - 1);
			}
			void swap(chained_hash_set &other) {
				std::swap(m_sizes_index,other.m_sizes_index);
				std::swap(m_length,other.m_length);
				m_container.swap(other.m_container);
			}
			size_t size() const {
				return m_length;
			}
			// Increase size to at least new_size buckets.
			void resize(const size_t &new_size) {
				const uint8 new_sizes_index = find_upper_size_index(new_size);
				if (new_sizes_index <= m_sizes_index) {
					__PDEBUG(std::cout << "Invalid resize requested, doing nothing.\n");
					return;
				}
				change_size(new_sizes_index);
			}
			void erase(const iterator &it) {
				p_assert(it.m_ht == this);
				bucket_type &bucket = m_container[it.m_vector_index];
				const size_t bucket_size = bucket.size();
				p_assert(bucket_size > 0);
				p_assert(it.m_bucket_index < bucket_size);
				if (it.m_bucket_index + 1 != bucket_size) {
					// If the element is not at the last slot of the bucket,
					// swap it with the last element.
					swapper()(bucket[it.m_bucket_index], bucket[bucket_size - 1]);
				}
				// Erase the last element of the bucket.
				bucket.pop_back();
				--m_length;
				if (static_cast<double>(m_length) * shrink_trigger <
					settings::load_factor() * sizes[m_sizes_index] && m_sizes_index > min_size_index) {
					change_size(m_sizes_index - 1);
				}
			}
			bool empty() const {
				return (m_length == 0);
			}
			void clear() {
				chained_hash_set new_ht;
				new_ht.m_sizes_index = min_size_index;
				new_ht.m_container.resize(sizes[new_ht.m_sizes_index]);
				swap(new_ht);
			}
		private:
			static size_t get_position(const size_t &h, const size_t &size) {
				return (h % size);
			}
			static uint8 find_upper_size_index(const size_t &size) {
				for (uint8 retval = 0; retval < sizes_size; ++retval) {
					if (sizes[retval] >= size) {
						return std::max<uint8>(min_size_index,retval);
					}
				}
				return min_size_index;
			}
			void change_size(const uint8 &new_index) {
				p_assert(new_index >= min_size_index);
				chained_hash_set new_ht;
				new_ht.m_sizes_index = new_index;
				new_ht.m_container.resize(sizes[new_ht.m_sizes_index]);
				const iterator it_f = end();
				for (iterator it = begin(); it != it_f; ++it) {
					new_ht.insert(*it);
				}
				swap(new_ht);
				__PDEBUG(std::cout << "Chained hash set changed to: " << sizes[m_sizes_index] << '\n');
			}
			// Increase the size of the container to the next size.
			void grow() {
				p_assert(m_sizes_index < sizes_size - 1);
				change_size(m_sizes_index + 1);
			}
		private:
			uint8				m_sizes_index;
			size_t				m_length;
			vector_type			m_container;
	};

	template <class T, class HashFunction, class EqualKey, class Allocator, class SwapFunction>
	const size_t chained_hash_set<T,HashFunction,EqualKey,Allocator,SwapFunction>::sizes[] = {
		1,
		3,
		5,
		11,
		23,
		53,
		97,
		193,
		389,
		769,
		1543,
		3079,
		6151,
		12289,
		24593,
		49157,
		98317,
		196613,
		393241,
		786433,
		1572869,
		3145739,
		6291469,
		12582917,
		25165843,
		50331653,
		100663319,
		201326611,
		402653189,
		805306457,
		1610612741,
		3221225473
#ifdef _PIRANHA_64BIT
		,
		6442450939,
		12884901893,
		25769803799,
		51539607551,
		103079215111,
		206158430209,
		412316860441,
		824633720831
#endif
	};
}

#endif
