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

#include <algorithm> // For std::max.
#include <exception> // For std::bad_alloc.
#include <vector>

#include "../config.h"
#include "../exceptions.h" // For out_of_memory.
#include "../math.h"
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
					void swap(term_type_ &t) {
						m_cf.swap(t.m_cf);
						int_swap(m_ckey,t.m_ckey);
					}
					mutable Cf	m_cf;
					Ckey		m_ckey;
			};
			template <int N>
			class bucket_type_ {
					p_static_check(N > 0 && (N == 1 || lg<N>::value > 0), "Invalid bucket size in cuckoo hash.");
				public:
					static const size_t size = static_cast<size_t>(N);
					bucket_type_() {
						for (size_t i = 0; i < size; ++i) {
							f[i] = false;
						}
					}
					term_type_	t[size];
					bool 		f[size];
			};
			// Bucket size = 4.
			typedef bucket_type_<4> bucket_type;
			static const size_t bsize = bucket_type::size;
			typedef typename Allocator::template rebind<bucket_type>::other allocator_type;
			typedef std::vector<bucket_type,allocator_type> container_type;
			static const size_t mults[4];
			static const size_t mults_size = 4;
#ifdef _PIRANHA_64BIT
			static const size_t sizes[64];
			static const size_t sizes_size = 64;
#else
			static const size_t sizes[32];
			static const size_t sizes_size = 32;
#endif
			p_static_check(mults_size % 2 == 0, "Mults size must be a multiple of 2.");
		public:
			typedef term_type_ term_type;
			class iterator
			{
					friend class coded_series_cuckoo_hash_table;
					explicit iterator(const coded_series_cuckoo_hash_table *ht): m_ht(ht), m_vindex(0), m_bindex(0) {
						// Go to the first occupied term if the table is not empty and the first
						// term is not occupied.
						if (m_ht->m_container.size() > 0 && !m_ht->m_container[0].f[0]) {
							next();
						}
					}
					explicit iterator(const coded_series_cuckoo_hash_table *ht, const size_t &vpos, const size_t &bpos):
						m_ht(ht), m_vindex(vpos), m_bindex(bpos) {
						p_assert(m_bindex < bsize);
					}
				public:
					iterator &operator++() {
						next();
						return *this;
					}
					const term_type &operator*() const {
						p_assert(m_vindex < m_ht->m_container.size());
						p_assert(m_bindex < bsize);
						p_assert(m_ht->m_container[m_vindex].f[m_bindex]);
						return m_ht->m_container[m_vindex].t[m_bindex];
					}
					const term_type *operator->() const {
						p_assert(m_vindex < m_ht->m_container.size());
						p_assert(m_bindex < bsize);
						p_assert(m_ht->m_container[m_vindex].f[m_bindex]);
						return &m_ht->m_container[m_vindex].t[m_bindex];
					}
					bool operator==(const iterator &it2) const {
						return (m_ht == it2.m_ht && m_vindex == it2.m_vindex && m_bindex == it2.m_bindex);
					}
					bool operator!=(const iterator &it2) const {
						return !(*this == it2);
					}
				private:
					void next() {
						const size_t vector_size = m_ht->m_container.size();
						size_t tmp_vindex = m_vindex, tmp_bindex = m_bindex;
						do {
							if (tmp_bindex == bsize - 1) {
								tmp_bindex = 0;
								++tmp_vindex;
							} else {
								++tmp_bindex;
							}
						} while (tmp_vindex < vector_size && !m_ht->m_container[tmp_vindex].f[tmp_bindex]);
						m_vindex = tmp_vindex;
						m_bindex = tmp_bindex;
					}
				private:
					const coded_series_cuckoo_hash_table	*m_ht;
					size_t									m_vindex;
					size_t									m_bindex;
			};
			coded_series_cuckoo_hash_table(): m_sizes_index(0), m_mults_index(0), m_length(0) {
				m_container.resize(sizes[2]);
				m_sizes_index = 2;
			}
			coded_series_cuckoo_hash_table(const size_t &size): m_mults_index(0), m_length(0) {
				uint8 sizes_index = find_upper_pow2_index(size / bsize);
				p_assert(sizes_index >= 2);
				while (true) {
					try {
						m_container.resize(sizes[sizes_index]);
						m_sizes_index = sizes_index;
						break;
					} catch (const out_of_memory &) {
						--sizes_index;
						if (sizes_index < 2) {
							throw out_of_memory("Not enough available memory to allocate a cuckoo hash table of size 4.");
						}
					} catch (const std::bad_alloc &) {
						--sizes_index;
						if (sizes_index < 2) {
							throw out_of_memory("Not enough physical memory to allocate a cuckoo hash table of size 4.");
						}
					}
				}
			}
			~coded_series_cuckoo_hash_table() {
				__PDEBUG(
				size_t i = 0;
				for (iterator it = begin(); it != end(); ++it) {
					++i;
				}
				p_assert(i == size());
				std::cout << "No problems. Final vector size: "<< m_container.size() << '\n';
				)
				p_assert(sizes[m_sizes_index] == m_container.size());
			}
			iterator begin() const {
				p_assert(sizes[m_sizes_index] == m_container.size());
				return iterator(this);
			}
			iterator end() const {
				p_assert(sizes[m_sizes_index] == m_container.size());
				return iterator(this, m_container.size(), 0);
			}
			iterator find(const Ckey &ckey) const {
				p_assert(sizes[m_sizes_index] == m_container.size());
				const size_t h = static_cast<size_t>(ckey), pos1 = position1(h), pos2 = position2(h);
				// TODO: replace with bit twiddling to reduce branching?
				for (size_t i = 0; i < bsize; ++i) {
					if (m_container[pos1].f[i] && m_container[pos1].t[i].m_ckey == ckey) {
						return iterator(this,pos1,i);
					}
				}
				for (size_t i = 0; i < bsize; ++i) {
					if (m_container[pos2].f[i] && m_container[pos2].t[i].m_ckey == ckey) {
						return iterator(this,pos2,i);
					}
				}
				return end();
			}
			size_t size() const {
				p_assert(sizes[m_sizes_index] == m_container.size());
				return m_length;
			}
			void insert(const term_type &t) {
				if ((static_cast<double>(m_length) + 1) >=
					settings::load_factor() * (sizes[m_sizes_index] * bsize)) {
					__PDEBUG(std::cout << "Max load factor exceeded, resizing." << '\n');
					increase_size();
				}
				term_type tmp_term;
				if (!attempt_insertion(t,tmp_term)) {
					if (!rehash()) {
						// If rehash was not successful, resize the table.
						increase_size();
					}
					// We still have to insert the displaced term that was left out from the failed attempt.
					// This is stored in the temporary term. We have a recursion going on here, it should
					// not matter much because most likely it is overpowered by the resize above. Probably
					// it could be turned into an iteration with some effort?
					insert(tmp_term);
				}
				p_assert(sizes[m_sizes_index] == m_container.size());
			}
			void swap(coded_series_cuckoo_hash_table &other) {
				int_swap(m_mults_index,other.m_mults_index);
				int_swap(m_sizes_index,other.m_sizes_index);
				int_swap(m_length,other.m_length);
				m_container.swap(other.m_container);
				p_assert(sizes[m_sizes_index] == m_container.size());
			}
		private:
			static uint8 find_upper_pow2_index(const size_t &size) {
				for (uint8 retval = 0; retval < sizes_size; ++retval) {
					if (sizes[retval] >= size) {
						return std::max<uint8>(static_cast<uint8>(2),retval);
					}
				}
				return static_cast<uint8>(2);
			}
			void increase_size() {
				coded_series_cuckoo_hash_table new_ht;
				new_ht.m_container.resize(sizes[m_sizes_index + 1]);
				new_ht.m_sizes_index = m_sizes_index + 1;
				term_type tmp_term;
				iterator it = begin();
				const iterator it_f = end();
				while (it != it_f) {
					if (!new_ht.attempt_insertion(*it,tmp_term)) {
						// TODO: here we should check with other hash functions before giving up and increasing size.
						__PDEBUG(std::cout << "Cuckoo hash table resize triggered during resize." << '\n');
						const size_t new_index = new_ht.m_sizes_index + 1;
						__PDEBUG(std::cout << "Next size: " << sizes[new_index] << '\n');
						// We do things in this order because if resizing fails due to insufficient memory,
						// then we must be sure that the size index is consistent with the initial size (otherwise
						// assertions can fail in the dtor, for example, or elsewhere).
						new_ht.m_container.clear();
						new_ht.m_sizes_index = 0;
						new_ht.m_container.resize(sizes[new_index]);
						new_ht.m_sizes_index = new_index;
						new_ht.m_length = 0;
						it = begin();
					} else {
						++it;
					}
				}
				swap(new_ht);
				p_assert(sizes[m_sizes_index] == m_container.size());
			}
			bool rehash() {
				for (uint8 new_mults_index = m_mults_index + 2; new_mults_index < mults_size; new_mults_index += 2) {
					coded_series_cuckoo_hash_table new_ht;
					new_ht.m_mults_index = new_mults_index;
					new_ht.m_container.resize(sizes[m_sizes_index]);
					new_ht.m_sizes_index = m_sizes_index;
					term_type tmp_term;
					const iterator it_f = end();
					for (iterator it = begin(); it != it_f; ++it) {
						if (!new_ht.attempt_insertion(*it,tmp_term)) {
							__PDEBUG(std::cout << "Cuckoo hash table insertion failure during rehash, "
								"will try next mult.\n");
							break;
						}
					}
					if (new_ht.size() == size()) {
						__PDEBUG(std::cout << "Rehash successful after " << ((new_mults_index - m_mults_index)/2) <<
							" tries\n";)
						// This means that we were able to insert all terms. Swap and return true.
						swap(new_ht);
						p_assert(sizes[m_sizes_index] == m_container.size());
						return true;
					}
					__PDEBUG(std::cout << "Rehash stopped at " << new_ht.size() << " out of " << size() << ".\n");
				}
				__PDEBUG(std::cout << "Mults exhausted, rehash failed.\n")
				p_assert(sizes[m_sizes_index] == m_container.size());
				return false;
			}
			bool attempt_insertion(const term_type &t, term_type &tmp_term) {
				const size_t h = static_cast<size_t>(t.m_ckey);
				size_t pos = position1(h);
				for (size_t i = 0; i < bsize; ++i) {
					// There's space in the bucket, rejoice!
					if (!m_container[pos].f[i]) {
						m_container[pos].f[i] = true;
						m_container[pos].t[i] = t;
						++m_length;
						p_assert(sizes[m_sizes_index] == m_container.size());
						return true;
					}
				}
				// No space was found in the first-choice bucket. Choose randomly(?) the index of the element
				// in the bucket that will be displaced.
				size_t dindex;
				if (bsize == 1) {
					// If we have single-term bucket, result will always be 0.
					dindex = 0;
				} else {
					// We know bucket size is power of 2 when it is not 1, use bit twiddling for mod.
					dindex = (pos & (bsize - 1));
				}
				tmp_term = m_container[pos].t[dindex];
				m_container[pos].t[dindex] = t;
				size_t counter = 0;
				while (swap_and_displace(tmp_term,pos)) {
					++counter;
					if (counter > 10) {
						__PDEBUG(std::cout << "Cuckoo loop detected, returning false. Load factor is: " <<
							((static_cast<double>(m_length) + 1) / (sizes[m_sizes_index] * bsize)) << '\n');
						p_assert(sizes[m_sizes_index] == m_container.size());
						return false;
					}
				}
				p_assert(sizes[m_sizes_index] == m_container.size());
				return true;
			}
			// Place tmp_term into its location other than orig_location, displacing, if necessary. an existing
			// term. If displacement takes place, retval will be true and the content of tmp_term
			// will be the displaced one. Otherwise return false.
			bool swap_and_displace(term_type &tmp_term, size_t &orig_location) {
				//__PDEBUG(std::cout << "Performing swap & displace." << '\n');
				// First thing we need to know if the original location was given by hash1 or hash2.
				const size_t h = static_cast<size_t>(tmp_term.m_ckey), pos1 = position1(h), pos2 = position2(h);
				size_t new_pos;
				if (orig_location == pos1) {
					// Original location was pos1, we want to move to pos2.
					new_pos = pos2;
				} else {
					// Original location was pos2, we want to move to pos1.
					new_pos = pos1;
				}
				for (size_t i = 0; i < bsize; ++i) {
					if (!m_container[new_pos].f[i]) {
						// Place found, rejoice!
						m_container[new_pos].t[i].swap(tmp_term);
						m_container[new_pos].f[i] = true;
						++m_length;
						p_assert(sizes[m_sizes_index] == m_container.size());
						return false;
					}
				}
				// No space was found in the first-choice bucket. Choose randomly(?) the index of the element
				// in the bucket that will be displaced.
				size_t dindex;
				if (bsize == 1) {
					// If we have single-term bucket, result will always be 0.
					dindex = 0;
				} else {
					// We know bucket size is power of 2 when it is not 1, use bit twiddling for mod.
					dindex = (new_pos & (bsize - 1));
				}
				// Displace the selected term.
				m_container[new_pos].t[dindex].swap(tmp_term);
				orig_location = new_pos;
				p_assert(sizes[m_sizes_index] == m_container.size());
				return true;
			}
			static size_t m_hash(const size_t &h, const size_t &mult, const size_t &sizes_index) {
				p_assert(mult * h < sizes[sizes_index]);
				return (mult * h) >> ((sizeof(max_fast_int) << 3) - sizes_index);
			}
			size_t position1(const size_t &h) const {
				return m_hash(h,mults[m_mults_index],m_sizes_index);
			}
			size_t position2(const size_t &h) const {
				return m_hash(h,mults[m_mults_index + 1],m_sizes_index);
			}
		private:
			uint8				m_sizes_index;
			uint8				m_mults_index;
			size_t				m_length;
			container_type		m_container;
	};

	template <class Cf, class Ckey, class Allocator>
	const size_t coded_series_cuckoo_hash_table<Cf,Ckey,Allocator>::sizes[] = {
		1,
		2,
		4,
		8,
		16,
		32,
		64,
		128,
		256,
		512,
		1024,
		2048,
		4096,
		8192,
		16384,
		32768,
		65536,
		131072,
		262144,
		524288,
		1048576,
		2097152,
		4194304,
		8388608,
		16777216,
		33554432,
		67108864,
		134217728,
		268435456,
		536870912,
		1073741824,
		2147483648u
#ifdef _PIRANHA_64BIT
		,
		4294967296,
		8589934592,
		17179869184,
		34359738368,
		68719476736,
		137438953472,
		274877906944,
		549755813888,
		1099511627776,
		2199023255552,
		4398046511104,
		8796093022208,
		17592186044416,
		35184372088832,
		70368744177664,
		140737488355328,
		281474976710656,
		562949953421312,
		1125899906842624,
		2251799813685248,
		4503599627370496,
		9007199254740992,
		18014398509481984,
		36028797018963968,
		72057594037927936,
		144115188075855872,
		288230376151711744,
		576460752303423488,
		1152921504606846976,
		2305843009213693952,
		4611686018427387904,
		9223372036854775808u
#endif
	};

#ifdef _PIRANHA_64BIT
#define MAX (18446744073709551616.)
#else
#define MAX (4294967296.)
#endif

	template <class Cf, class Ckey, class Allocator>
	const size_t coded_series_cuckoo_hash_table<Cf,Ckey,Allocator>::mults[] = {
		static_cast<size_t>(.7320508075688772 * MAX),
		static_cast<size_t>(.2360679774997898 * MAX),
// 		.6180339887498949 * MAX,
// 		.4658204617032757 * MAX,
		static_cast<size_t>(.6457513110645907 * MAX),
		static_cast<size_t>(.3166247903553998 * MAX)/*,
		.6055512754639891 * MAX,
		.1231056256176606 * MAX,
		.3588989435406740 * MAX,
		.7958315233127191 * MAX,
		.3851648071345037 * MAX,
		.5677643628300215 * MAX*/
	};
#undef MAX
}

#endif
