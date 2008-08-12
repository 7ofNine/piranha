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

#ifndef PIRANHA_CODED_SERIES_CUCKOO_HASH_TABLE_H
#define PIRANHA_CODED_SERIES_CUCKOO_HASH_TABLE_H

#include <algorithm> // For swap.
#include <boost/array.hpp>
#include <boost/functional/hash.hpp>
#include <vector>

#include "../config.h"
#include "../memory.h"
#include "../p_assert.h"
#include "../settings.h"

namespace piranha
{
	template <class Cf, class Ckey, class Allocator>
	class coded_series_cuckoo_hash_table
	{
			class term_type_ {
				public:
					term_type_():m_cf(),m_ckey() {}
					term_type_(const Cf &cf, const Ckey &ckey):m_cf(cf),m_ckey(ckey) {}
					mutable Cf	m_cf;
					Ckey		m_ckey;
			};
			class bucket_type {
				public:
					term_type_	m_term;
					size_t		m_hash;
					bool		m_flag;
			};
			typedef typename Allocator::template rebind<bucket_type>::other allocator_type;
			typedef std::vector<bucket_type,allocator_type> container_type;
			typedef boost::array<size_t,37> sizes_vector_type;
			static const sizes_vector_type sizes;
		public:
			typedef term_type_ term_type;
			class iterator
			{
					friend class coded_series_cuckoo_hash_table;
					iterator(const coded_series_cuckoo_hash_table *ht): m_ht(ht), m_index(0) {
						// Go to the first occupied bucket if the table is not empty and the first
						// bucket is not occupied.
						if (m_ht->m_container.size() > 0 && !m_ht->m_container[0].m_flag) {
							next();
						}
					}
					iterator(const coded_series_cuckoo_hash_table *ht, const size_t &pos): m_ht(ht), m_index(pos) {}
				public:
					iterator &operator++() {
						next();
						return *this;
					}
					const term_type &operator*() const {
						p_assert(m_index < m_ht->m_container.size());
						p_assert(m_ht->m_container[m_index].m_flag);
						return m_ht->m_container[m_index].m_term;
					}
					const term_type *operator->() const {
						p_assert(m_index < m_ht->m_container.size());
						p_assert(m_ht->m_container[m_index].m_flag);
						return &m_ht->m_container[m_index].m_term;
					}
					bool operator==(const iterator &it2) const {
						return (m_ht == it2.m_ht && m_index == it2.m_index);
					}
					bool operator!=(const iterator &it2) const {
						return !(*this == it2);
					}
				private:
					void next() {
						const size_t vector_size = sizes[m_ht->m_sizes_index];
						size_t tmp_index = m_index + 1;
						while (tmp_index < vector_size && !m_ht->m_container[tmp_index].m_flag) {
							++tmp_index;
						}
						m_index = tmp_index;
					}
				private:
					const coded_series_cuckoo_hash_table	*m_ht;
					size_t									m_index;
			};
			coded_series_cuckoo_hash_table(): m_length(0), m_sizes_index(0), m_container(sizes[0]) {}
			~coded_series_cuckoo_hash_table() {
				__PDEBUG(std::cout << "On destruction, the vector size of coded_series_cuckoo_hash_table was: "
								   << m_container.size() << '\n');
			}
			iterator begin() const {
				return iterator(this);
			}
			iterator end() const {
				p_assert(sizes[m_sizes_index] == m_container.size());
				return iterator(this, sizes[m_sizes_index]);
			}
			iterator find(const Ckey &ckey) const {
				p_assert(sizes[m_sizes_index] == m_container.size());
				const size_t v_size = sizes[m_sizes_index], h1 = static_cast<size_t>(ckey), pos1 = h1 % v_size;
				// TODO: replace with bit twiddling to reduce branching?
				if (m_container[pos1].m_flag && m_container[pos1].m_term.m_ckey == ckey) {
					return iterator(this,pos1);
				}
				const size_t pos2 = (~h1) % v_size;
				if (m_container[pos2].m_flag && m_container[pos2].m_term.m_ckey == ckey) {
					return iterator(this,pos2);
				}
				return end();
			}
			size_t size() const {
				return m_length;
			}
			void insert(const term_type &t) {
				if (((m_length + 1) << 1) >= sizes[m_sizes_index]) {
					__PDEBUG(std::cout << "Load factor exceeded, resizing." << '\n');
					increase_size();
				}
				bucket_type tmp_bucket;
				if (!attempt_insertion(t,tmp_bucket)) {
					// If we fail insertion, we must increase size.
					increase_size();
					// We still have to insert the displaced term that was left out from the failed attempt.
					// This is stored in the temporary bucket. We have a recursion going on here, it should
					// not matter much because most likely it is overpowered by the resize above. Probably
					// it could be turned into an iteration with some effort?
					insert(tmp_bucket.m_term);
				}
			}
		private:
			void increase_size() {
				coded_series_cuckoo_hash_table new_ht;
				new_ht.m_sizes_index = m_sizes_index + 1;
				new_ht.m_container.resize(sizes[new_ht.m_sizes_index]);
				bucket_type tmp_bucket;
				iterator it = begin();
				const iterator it_f = end();
				while (it != it_f) {
					if (!new_ht.attempt_insertion(*it,tmp_bucket)) {
						__PDEBUG(std::cout << "Cuckoo hash table resize triggered during resize." << '\n');
						++new_ht.m_sizes_index;
						__PDEBUG(std::cout << "New size: " << sizes[new_ht.m_sizes_index] << '\n');
						new_ht.m_container.clear();
						new_ht.m_container.resize(sizes[new_ht.m_sizes_index]);
						new_ht.m_length = 0;
						it = begin();
					} else {
						++it;
					}
				}
				m_container.swap(new_ht.m_container);
				m_length = new_ht.m_length;
				m_sizes_index = new_ht.m_sizes_index;
// std::cout << "After resize:\n";
// for (iterator it = begin(); it != end(); ++it) {
// 	std::cout << it->m_ckey << '\n';
// }
// std::cout << "----------\n";
			}
			bool attempt_insertion(const term_type &t, bucket_type &tmp_bucket) {
				const size_t vector_size = sizes[m_sizes_index], h = static_cast<size_t>(t.m_ckey),
					pos = h % vector_size;
				if (m_container[pos].m_flag) {
					// Current bucket is kicked out into tmp_bucket, t takes its place
					tmp_bucket = m_container[pos];
					m_container[pos].m_term = t;
					m_container[pos].m_hash = h;
					size_t counter = 0;
					while (swap_and_displace(tmp_bucket,vector_size)) {
						++counter;
						if (counter > 19) {
// std::cout << "tmp contains: " << tmp_bucket.m_term.m_ckey << '\n';
// std::cout << "size is: " << size() << '\n';
// for (iterator it = begin(); it != end(); ++it) {
// 	std::cout << it->m_ckey << '\n';
// }
							__PDEBUG(std::cout << "Cuckoo loop detected, will increase size and rebuild.\n");
							return false;
						}
					}
				} else {
					m_container[pos].m_flag = true;
					m_container[pos].m_term = t;
					m_container[pos].m_hash = h;
				}
				++m_length;
				return true;
			}
			// Place tmp_bucket into its alternative location, displacing, if necessary. an existing
			// bucket. If displacement takes place, retval will be true and the content of tmp_bucket
			// will be the displaced one. Otherwise return false.
			bool swap_and_displace(bucket_type &tmp_bucket, const size_t &vector_size) {
				//__PDEBUG(std::cout << "Performing swap & displace." << '\n');
				const size_t alt_h = ~tmp_bucket.m_hash, alt_pos = alt_h % vector_size;
				tmp_bucket.m_hash = alt_h;
				const bool retval = m_container[alt_pos].m_flag;
				if (retval) {
#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
					m_container[alt_pos].m_term.m_cf.swap(tmp_bucket.m_term.m_cf);
					SWAP(m_container[alt_pos].m_term.m_ckey,tmp_bucket.m_term.m_ckey);
					SWAP(m_container[alt_pos].m_hash,tmp_bucket.m_hash);
					SWAP(m_container[alt_pos].m_flag,tmp_bucket.m_flag);
#undef SWAP
				} else {
					m_container[alt_pos] = tmp_bucket;
				}
				p_assert(m_container[alt_pos].m_hash % vector_size == alt_pos);
				return retval;
			}
		private:
			size_t				m_length;
			size_t				m_sizes_index;
			container_type		m_container;
	};

	template <class Cf, class Ckey, class Allocator>
	const typename coded_series_cuckoo_hash_table<Cf,Ckey,Allocator>::sizes_vector_type
		coded_series_cuckoo_hash_table<Cf,Ckey,Allocator>::sizes = { {
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
	} };
}

#endif
