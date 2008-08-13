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

#include <boost/array.hpp>
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
			typedef typename Allocator::template rebind<term_type_>::other allocator_type;
			typedef std::vector<term_type_,allocator_type> container_type;
			typedef boost::array<size_t,37> sizes_vector_type;
			static const sizes_vector_type sizes;
		public:
			typedef term_type_ term_type;
			class iterator
			{
					friend class coded_series_cuckoo_hash_table;
					iterator(const coded_series_cuckoo_hash_table *ht): m_ht(ht), m_index(0) {
						// Go to the first occupied term if the table is not empty and the first
						// term is not occupied.
						if (m_ht->m_container.size() > 0 && !m_ht->m_flags[0]) {
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
						p_assert(m_ht->m_flags[m_index]);
						return m_ht->m_container[m_index];
					}
					const term_type *operator->() const {
						p_assert(m_index < m_ht->m_container.size());
						p_assert(m_ht->m_flags[m_index]);
						return &m_ht->m_container[m_index];
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
						while (tmp_index < vector_size && !m_ht->m_flags[tmp_index]) {
							++tmp_index;
						}
						m_index = tmp_index;
					}
				private:
					const coded_series_cuckoo_hash_table	*m_ht;
					size_t									m_index;
			};
			coded_series_cuckoo_hash_table(): m_length(0), m_sizes_index(0), m_container(sizes[0]), m_flags(sizes[0]) {}
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
				const size_t v_size = sizes[m_sizes_index], h1 = hash1(ckey), pos1 = h1 % v_size;
				// TODO: replace with bit twiddling to reduce branching?
				if (m_flags[pos1] && m_container[pos1].m_ckey == ckey) {
					return iterator(this,pos1);
				}
				const size_t pos2 = hash2(ckey) % v_size;
				if (m_flags[pos2] && m_container[pos2].m_ckey == ckey) {
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
				term_type tmp_term;
				if (!attempt_insertion(t,tmp_term)) {
					// If we fail insertion, we must increase size.
					increase_size();
					// We still have to insert the displaced term that was left out from the failed attempt.
					// This is stored in the temporary term. We have a recursion going on here, it should
					// not matter much because most likely it is overpowered by the resize above. Probably
					// it could be turned into an iteration with some effort?
					insert(tmp_term);
				}
			}
		private:
			void increase_size() {
				coded_series_cuckoo_hash_table new_ht;
				new_ht.m_sizes_index = m_sizes_index + 1;
				new_ht.m_container.resize(sizes[new_ht.m_sizes_index]);
				new_ht.m_flags.resize(sizes[new_ht.m_sizes_index]);
				term_type tmp_term;
				iterator it = begin();
				const iterator it_f = end();
				while (it != it_f) {
					if (!new_ht.attempt_insertion(*it,tmp_term)) {
						__PDEBUG(std::cout << "Cuckoo hash table resize triggered during resize." << '\n');
						++new_ht.m_sizes_index;
						__PDEBUG(std::cout << "New size: " << sizes[new_ht.m_sizes_index] << '\n');
						new_ht.m_container.clear();
						new_ht.m_container.resize(sizes[new_ht.m_sizes_index]);
						new_ht.m_flags.clear();
						new_ht.m_flags.resize(sizes[new_ht.m_sizes_index]);
						new_ht.m_length = 0;
						it = begin();
					} else {
						++it;
					}
				}
				m_container.swap(new_ht.m_container);
				m_flags.swap(new_ht.m_flags);
				m_length = new_ht.m_length;
				m_sizes_index = new_ht.m_sizes_index;
			}
			bool attempt_insertion(const term_type &t, term_type &tmp_term) {
				const size_t vector_size = sizes[m_sizes_index], h = static_cast<size_t>(t.m_ckey),
					pos = h % vector_size;
				if (m_flags[pos]) {
					// Current term is kicked out into tmp_term, t takes its place
					tmp_term = m_container[pos];
					m_container[pos] = t;
					size_t counter = 0;
					while (swap_and_displace(tmp_term,pos,vector_size)) {
						++counter;
						if (counter > 19) {
							__PDEBUG(std::cout << "Cuckoo loop detected, will increase size and rebuild.\n");
							return false;
						}
					}
				} else {
					m_flags[pos] = true;
					m_container[pos] = t;
				}
				++m_length;
				return true;
			}
			// Place tmp_term into its location other than orig_location, displacing, if necessary. an existing
			// term. If displacement takes place, retval will be true and the content of tmp_term
			// will be the displaced one. Otherwise return false.
			bool swap_and_displace(term_type &tmp_term, const size_t &orig_location, const size_t &vector_size) {
				//__PDEBUG(std::cout << "Performing swap & displace." << '\n');
				const size_t pos1 = hash1(tmp_term.m_ckey) % vector_size;
				size_t new_pos;
				if (orig_location == pos1) {
					// Original location was pos1, we want to move to pos2.
					new_pos = hash2(tmp_term.m_ckey) % vector_size;;
				} else {
					// Original location was pos2, we want to move to pos1.
					new_pos = pos1;
				}
				if (m_flags[new_pos]) {
					// We have to displace an existing term.
#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
					m_container[new_pos].m_cf.swap(tmp_term.m_cf);
					SWAP(m_container[new_pos].m_ckey,tmp_term.m_ckey);
					return true;
#undef SWAP
				} else {
					// Destination is not taken, occupy it.
					m_container[new_pos] = tmp_term;
					m_flags[new_pos] = true;
					return false;
				}
			}
			static size_t hash1(const Ckey &ckey) {
				//return static_cast<size_t>(ckey & __MASK);
				return static_cast<size_t>(ckey);
			}
			static size_t hash2(const Ckey &ckey) {
				//return static_cast<size_t>((ckey >> 32) & __MASK);
				//return (static_cast<size_t>(ckey) + 0x9e3779b9 + (ckey<<6) + (ckey>>2));
				return ~hash1(ckey);
			}
		private:
			size_t				m_length;
			size_t				m_sizes_index;
			container_type		m_container;
			std::vector<char>	m_flags;
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
