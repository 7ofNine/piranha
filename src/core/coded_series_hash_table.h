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

#include <algorithm>
#include <utility> // For std::pair.

#include "config.h"
#include "p_assert.h"
#include "settings.h"

namespace piranha
{
	template <class T, class Allocator>
	class coded_series_hash_table
	{
			template <int N>
			class bucket_type_
			{
				public:
					bucket_type_() {
						for (size_t i = 0; i < N; ++i) {
							f[i] = false;
						}
					}
					T		t[N];
					bool	f[N];
			};
			// Configuration options.
			static const size_t bucket_size		= 4;
			static const size_t min_size_index	= 0;
			// Number of extra buckets.
			static const size_t neb				= 128;
			// Configuration options end here.
			static const size_t sizes_size =
#ifdef _PIRANHA_64BIT
				40;
#else
				32;
#endif
			static const size_t sizes[sizes_size];
			typedef bucket_type_<bucket_size> bucket_type;;
			typedef bucket_type *container_type;
			p_static_check(bucket_size > 0, "Size of bucket is not strictly positive.");
			p_static_check(min_size_index >= 0, "min_size_index must be non-negative.");
		public:
			typedef T key_type;
			typedef counting_allocator<bucket_type,Allocator> allocator_type;
			class iterator
			{
					friend class coded_series_hash_table;
					iterator(const coded_series_hash_table *p): m_ht(p), m_vector_index(0), m_bucket_index(0) {
						// If the first slot is not taken, find the next one.
						if (!m_ht->m_container[m_vector_index].f[m_bucket_index]) {
							next();
						}
					}
					iterator(const coded_series_hash_table *p, const size_t &vi, const size_t &bi):
							m_ht(p), m_vector_index(vi), m_bucket_index(bi) {}
				public:
					iterator &operator++() {
						next();
						return *this;
					}
					const key_type &operator*() const {
						p_assert(m_vector_index < sizes[m_ht->m_size_index] + neb);
						p_assert(m_bucket_index < bucket_size);
						p_assert(m_ht->m_container[m_vector_index].f[m_bucket_index]);
						return m_ht->m_container[m_vector_index].t[m_bucket_index];
					}
					const key_type *operator->() const {
						p_assert(m_vector_index < sizes[m_ht->m_size_index] + neb);
						p_assert(m_bucket_index < bucket_size);
						p_assert(m_ht->m_container[m_vector_index].f[m_bucket_index]);
						return &m_ht->m_container[m_vector_index].t[m_bucket_index];
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
						const size_t vector_size = sizes[m_ht->m_size_index] + neb;
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
							if (m_vector_index == vector_size || m_ht->m_container[m_vector_index].f[m_bucket_index]) {
								break;
							}
						}
					}
				private:
					const coded_series_hash_table	*m_ht;
					size_t							m_vector_index;
					size_t							m_bucket_index;
			};
			coded_series_hash_table(): m_size_index(min_size_index),m_length(0) {
				init();
			}
			coded_series_hash_table(const size_t &size_hint):
				m_size_index(find_upper_size_index(size_hint / bucket_size + 1)),m_length(0) {
				init();
			}
			~coded_series_hash_table() {
				__PDEBUG(std::cout << "On destruction, the vector size of coded_series_hash_table was: "
								   << sizes[m_size_index] << '\n');
				destroy();
			}
			iterator begin() const {
				return iterator(this);
			}
			iterator end() const {
				return iterator(this, sizes[m_size_index] + neb, 0);
			}
			std::pair<bool,iterator> find(const key_type &key) const {
				const size_t vector_size = sizes[m_size_index], vector_pos = key.hash_value() % vector_size;
				p_assert(vector_pos < sizes[m_size_index]);
				const bucket_type &bucket = m_container[vector_pos];
				// Now examine all elements in the bucket.
				for (size_t i = 0; i < bucket_size; ++i) {
					// If the slot in the bucket is not taken (which means there are no more elements to examine),
					// it means that key was not found.
					if (!bucket.f[i]) {
						return std::make_pair(false,iterator(this, vector_pos, i));
					}
					if (bucket.t[i] == key) {
						// If we found an occupied bucket slot, examine the key to see whether it matches or not
						// with t's.
						// If it does not match, let's move to the next bucket element.
						return std::make_pair(true,iterator(this, vector_pos, i));
					}
				}
				// Examined all the elements in the destination bucket. Examine the extra buckets.
				const bucket_type *extra_bucket = m_container + vector_size;
				for (size_t b = 0; b < neb; ++b) {
					for (size_t i = 0; i < bucket_size; ++i) {
						if (!extra_bucket[b].f[i]) {
							return std::make_pair(false,iterator(this, vector_size + b, i));
						}
						if (extra_bucket[b].t[i] == key) {
							return std::make_pair(true,iterator(this, vector_size + b, i));
						}
					}
				}
				// All the elements of the extra buckets were taken, we examined them but found no match.
				return std::make_pair(false,iterator(this, vector_size + neb, bucket_size));
			}
			void insert(const key_type &key, const iterator &it) {
				if (!attempt_insertion(key,it)) {
					iterator tmp = it;
					do {
						__PDEBUG(std::cout << "Started resizing coded series hash table." << '\n');
						increase_size();
						tmp = find(key).second;
						__PDEBUG(std::cout << "Resized coded series hash table." << '\n');
					} while (!attempt_insertion(key,tmp));
				}
			}
			size_t size() const {
				return m_length;
			}
		private:
			static size_t find_upper_size_index(const size_t &size) {
				for (size_t retval = 0; retval < sizes_size; ++retval) {
					if (sizes[retval] >= size) {
						return std::max<size_t>(min_size_index,retval);
					}
				}
				return min_size_index;
			}
			// Allocate and default-construct according to m_size_index.
			void init() {
				const size_t size = sizes[m_size_index] + neb;
				const bucket_type bucket;
				allocator_type a;
				m_container = a.allocate(size);
				for (size_t i = 0; i < size; ++i) {
					a.construct(m_container + i, bucket);
				}
			}
			// Destroy and deallocate.
			void destroy() {
				const size_t size = sizes[m_size_index] + neb;
				allocator_type a;
				for (size_t i = 0; i < size; ++i) {
					a.destroy(m_container + i);
				}
				a.deallocate(m_container, size);
			}
			bool attempt_insertion(const key_type &key, const iterator &it) {
				const size_t bucket_index = it.m_bucket_index;
				if (unlikely(bucket_index == bucket_size)) {
					return false;
				}
				p_assert(bucket_index < bucket_size);
				p_assert(it.m_vector_index < sizes[m_size_index] + neb);
				bucket_type &bucket = m_container[it.m_vector_index];
				p_assert(!bucket.f[bucket_index]);
				bucket.f[bucket_index] = true;
				bucket.t[bucket_index] = key;
				++m_length;
				return true;
			}
			// Increase size of the container to the next prime size.
			void increase_size() {
				__PDEBUG(std::cout << "Increase size requested at load factor: " <<
					((double)m_length / sizes[m_size_index]) << '\n');
				coded_series_hash_table new_ht;
				new_ht.destroy();
				new_ht.m_size_index = m_size_index + 2;
				new_ht.init();
				const iterator it_i = begin(), it_f = end();
				iterator it = it_i;
				while (it != it_f) {
					const std::pair<bool,iterator> res = new_ht.find(*it);
					p_assert(!res.first);
					if (!new_ht.attempt_insertion(*it,res.second)) {
						// NOTICE: here maybe we can use swapping instead of copying. The only problem is that
						// resizing can fail. In that case, we should swap back everything, if possible, and re-attempt
						// the resize with a bigger value.
						__PDEBUG(std::cout << "Hash table resize triggered during resize." << '\n');
						// TODO: add check for excessive size here. It must be here to avoid problems
						// with exception throwing.
						new_ht.destroy();
						++new_ht.m_size_index;
						new_ht.init();
						it = it_i;
					} else {
						++it;
					}
				}
				std::swap(m_container,new_ht.m_container);
				std::swap(m_size_index,new_ht.m_size_index);
				std::swap(m_length,new_ht.m_length);
			}
		private:
			size_t			m_size_index;
			size_t			m_length;
			container_type	m_container;
	};

	template <class T, class Allocator>
	const size_t coded_series_hash_table<T,Allocator>::min_size_index;

	template <class T, class Allocator>
	const size_t coded_series_hash_table<T,Allocator>::bucket_size;

	template <class T, class Allocator>
	const size_t coded_series_hash_table<T,Allocator>::sizes[] = {
		1u,
		3u,
		5u,
		11u,
		23u,
		53u,
		97u,
		193u,
		389u,
		769u,
		1543u,
		3079u,
		6151u,
		12289u,
		24593u,
		49157u,
		98317u,
		196613u,
		393241u,
		786433u,
		1572869u,
		3145739u,
		6291469u,
		12582917u,
		25165843u,
		50331653u,
		100663319u,
		201326611u,
		402653189u,
		805306457u,
		1610612741u,
		3221225473u
#ifdef _PIRANHA_64BIT
		,
		6442450939u,
		12884901893u,
		25769803799u,
		51539607551u,
		103079215111u,
		206158430209u,
		412316860441u,
		824633720831u
#endif
	};
}

#endif
