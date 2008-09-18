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

#ifndef PIRANHA_CUCKOO_HASH_SET_H
#define PIRANHA_CUCKOO_HASH_SET_H

#include <algorithm> // For std::max and std::swap.
#include <boost/functional/hash.hpp> // For default hash functor.
#include <boost/integer_traits.hpp> // For maximum values of size_t.
#include <exception> // For std::bad_alloc.
#include <functional>
#include <memory>
#include <utility> // For std::pair.
#include <sstream> // For building correct error message.
#include <string>

#include "common_functors.h"
#include "config.h" // For p_static_check.
#include "exceptions.h" // For out_of_memory.
#include "math.h"
#include "memory.h"
#include "p_assert.h"
#include "settings.h"

namespace piranha
{
	template <class T, class HashFunction = boost::hash<T>, class EqualKey = std::equal_to<T>,
	class Allocator = std::allocator<char>, class SwapFunction = std_swap<T> >
	class cuckoo_hash_set
	{
			template <int N>
			class bucket_type_ {
					p_static_check(N > 0 && (N == 1 || lg<N>::value > 0), "Invalid bucket size in cuckoo hash set.");
				public:
					static const size_t size = static_cast<size_t>(N);
					bucket_type_() {
						for (size_t i = 0; i < size; ++i) {
							f[i] = false;
						}
					}
					T		t[size];
					bool 	f[size];
			};
			// Configuration options.
			static const size_t bucket_size =			4;
			static const size_t mults_size = 			10;
			static const size_t max_rehash_tries =		1;
			static const size_t min_size_index =		1;
			static const size_t cuckoo_loop_threshold =	10;
			static const size_t shrink_trigger = 		10;
			// Configuration options stop here.
			static const size_t sizes_size =
#ifdef _PIRANHA_64BIT
				64;
#else
				32;
#endif
			typedef bucket_type_<bucket_size> bucket_type;
			static const size_t mults[mults_size];
			static const size_t sizes[sizes_size];
			p_static_check(mults_size % 2 == 0, "Mults size must be a multiple of 2.");
			p_static_check(shrink_trigger > 1, "Shrink trigger must be strictly greater than 1.");
			p_static_check(min_size_index > 0, "Min size index must be at least 1.");
		public:
			// Public typedefs.
			typedef T value_type;
			typedef T key_type;
			typedef HashFunction hasher;
			typedef EqualKey key_equal;
			typedef typename Allocator::template rebind<bucket_type>::other allocator_type;
			typedef SwapFunction swapper;
			typedef T * pointer;
			typedef const T * const_pointer;
			typedef T & reference;
			typedef const T & const_reference;
			typedef size_t size_type;
			typedef ptrdiff_t difference_type;
			// Iterator class.
			class iterator
			{
					friend class cuckoo_hash_set;
					explicit iterator(const cuckoo_hash_set *ht): m_ht(ht), m_vindex(0), m_bindex(0) {
						p_assert(sizes[m_ht->m_sizes_index] > 0);
						// Go to the first occupied slot if the first one isn't.
						if (!m_ht->m_container[0].f[0]) {
							next();
						}
					}
					explicit iterator(const cuckoo_hash_set *ht, const size_t &vpos, const size_t &bpos):
						m_ht(ht), m_vindex(vpos), m_bindex(bpos) {
						p_assert(m_bindex < bucket_size);
					}
				public:
					iterator &operator++() {
						next();
						return *this;
					}
					const_reference operator*() const {
						p_assert(m_vindex < sizes[m_ht->m_sizes_index]);
						p_assert(m_bindex < bucket_size);
						p_assert(m_ht->m_container[m_vindex].f[m_bindex]);
						return m_ht->m_container[m_vindex].t[m_bindex];
					}
					const_pointer operator->() const {
						p_assert(m_vindex < sizes[m_ht->m_sizes_index]);
						p_assert(m_bindex < bucket_size);
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
						const size_t vector_size = sizes[m_ht->m_sizes_index];
						const bucket_type *container = m_ht->m_container;
						size_t tmp_vindex = m_vindex, tmp_bindex = m_bindex;
						do {
							if (tmp_bindex == bucket_size - 1) {
								tmp_bindex = 0;
								++tmp_vindex;
							} else {
								++tmp_bindex;
							}
						} while (tmp_vindex < vector_size && !container[tmp_vindex].f[tmp_bindex]);
						m_vindex = tmp_vindex;
						m_bindex = tmp_bindex;
					}
				private:
					const cuckoo_hash_set	*m_ht;
					size_t					m_vindex;
					size_t					m_bindex;
			};
			typedef iterator const_iterator;
			// Constructors.
			cuckoo_hash_set(): m_sizes_index(min_size_index), m_mults_index(0), m_length(0) {
				init();
			}
			cuckoo_hash_set(const size_t &size): m_mults_index(0), m_length(0) {
				m_sizes_index = find_upper_size_index(size / bucket_size);
				p_assert(m_sizes_index >= min_size_index);
				while (true) {
					try {
						init();
						break;
					} catch (const out_of_memory &) {
						--m_sizes_index;
						if (m_sizes_index < min_size_index) {
							std::ostringstream stream;
							stream << "Not enough available memory to allocate a cuckoo hash set of size " <<
								sizes[min_size_index] << '.';
							throw out_of_memory(stream.str().c_str());
						}
					} catch (const std::bad_alloc &) {
						--m_sizes_index;
						if (m_sizes_index < min_size_index) {
							std::ostringstream stream;
							stream << "Not enough physical memory to allocate a cuckoo hash set of size " <<
								sizes[min_size_index] << '.';
							throw out_of_memory(stream.str().c_str());
						}
					}
				}
			}
			cuckoo_hash_set(const cuckoo_hash_set &c):
				m_sizes_index(c.m_sizes_index),
				m_mults_index(c.m_mults_index),
				m_length(c.m_length),
				m_container(allocator_type().allocate(sizes[m_sizes_index])) {
				// Memory has been allocated, now copy over from c.
				assign_from_other(c);
			}
			cuckoo_hash_set &operator=(const cuckoo_hash_set &c) {
				// Handle self-assignment.
				if (this == &c) {
					return *this;
				}
				// First let's auto-destroy.
				destroy();
				// Then assign the data members.
				m_sizes_index = c.m_sizes_index;
				m_mults_index = c.m_mults_index;
				m_length = c.m_length;
				// Finally, allocate the needed space and copy over from c.
				m_container = allocator_type().allocate(sizes[m_sizes_index]);
				assign_from_other(c);
				return *this;
			}
			~cuckoo_hash_set() {
				__PDEBUG(
				size_t i = 0;
				for (iterator it = begin(); it != end(); ++it) {
					++i;
				}
				p_assert(i == size());
				std::cout << "No problems. Final vector size: "<< sizes[m_sizes_index] << '\n';
				)
				destroy();
			}
			iterator begin() const {
				return iterator(this);
			}
			iterator end() const {
				return iterator(this, sizes[m_sizes_index], 0);
			}
			iterator find(const key_type &k) const {
				const size_t h = hasher()(k), pos1 = position1(h), pos2 = position2(h);
				bucket_type *container = m_container;
				// TODO: replace with bit twiddling to reduce branching?
				for (size_t i = 0; i < bucket_size; ++i) {
					if (container[pos1].f[i] && container[pos1].t[i] == k) {
						return iterator(this,pos1,i);
					}
				}
				for (size_t i = 0; i < bucket_size; ++i) {
					if (container[pos2].f[i] && container[pos2].t[i] == k) {
						return iterator(this,pos2,i);
					}
				}
				return end();
			}
			size_t size() const {
				return m_length;
			}
			bool empty() const {
				return (m_length == 0);
			}
			void insert(const key_type &k) {
				if ((static_cast<double>(m_length) + 1) >=
					settings::load_factor() * sizes[m_sizes_index] * bucket_size) {
					__PDEBUG(std::cout << "Max load factor exceeded, resizing." << '\n');
					grow();
				}
				key_type tmp_key;
				if (!attempt_insertion(k,tmp_key)) {
					if (!rehash()) {
						// If rehash was not successful, resize the table.
						grow();
					}
					// We still have to insert the displaced key that was left out from the failed attempt.
					// This is stored in the temporary key. We have a recursion going on here, it should
					// not matter much because most likely it is overpowered by the resize above. Probably
					// it could be turned into an iteration with some effort?
					insert(tmp_key);
				}
			}
			void erase(const iterator &it) {
				p_assert(it.m_ht == this);
				m_container[it.m_vindex].f[it.m_bindex] = false;
				m_container[it.m_vindex].t[it.m_bindex] = key_type();
				--m_length;
				if (static_cast<double>(m_length) * shrink_trigger < settings::load_factor() *
					sizes[m_sizes_index] * bucket_size) {
					attempt_shrink();
				}
			}
			void swap(cuckoo_hash_set &other) {
				std::swap(m_mults_index,other.m_mults_index);
				std::swap(m_sizes_index,other.m_sizes_index);
				std::swap(m_length,other.m_length);
				std::swap(m_container,other.m_container);
			}
		private:
			// Allocate and default-construct buckets according to the size specified by m_sizes_index.
			void init() {
				allocator_type a;
				m_container = a.allocate(sizes[m_sizes_index]);
				bucket_type *container = m_container;
				const bucket_type bucket;
				for (size_t i =  0; i < sizes[m_sizes_index]; ++i) {
					a.construct(container + i,bucket);
				}
			}
			// Destroy and deallocate buckets according to the size specified by m_sizes_index.
			void destroy() {
				allocator_type a;
				bucket_type *container = m_container;
				for (size_t i =  0; i < sizes[m_sizes_index]; ++i) {
					a.destroy(container + i);
				}
				a.deallocate(container,sizes[m_sizes_index]);
			}
			// Resize to new sizes index, destroying everything in the container.
			void resize(const uint8 &new_sizes_index) {
				destroy();
				m_sizes_index = new_sizes_index;
				m_length = 0;
				init();
			}
			// Construct buckets into m_container from other hash table c.
			void assign_from_other(const cuckoo_hash_set &c) {
				p_assert(m_sizes_index == c.m_sizes_index);
				allocator_type a;
				bucket_type *container = m_container;
				const bucket_type *container_other = c.m_container;
				const size_t size = sizes[m_sizes_index];
				for (size_t i = 0; i < size; ++i) {
					a.construct(container + i, container_other[i]);
				}
			}
			static uint8 find_upper_size_index(const size_t &size) {
				for (uint8 retval = 0; retval < sizes_size; ++retval) {
					if (sizes[retval] >= size) {
						return std::max<uint8>(min_size_index,retval);
					}
				}
				return min_size_index;
			}
			void grow() {
				if (size_t(m_sizes_index + 1) == sizes_size) {
					throw out_of_memory("Cuckoo hash table: growth exceeds the limits imposed for this architecture.");
				}
				cuckoo_hash_set new_ht;
				new_ht.resize(m_sizes_index + 1);
				// Assign current mults to the new table.
				new_ht.m_mults_index = m_mults_index;
				key_type tmp_key;
				iterator it = begin();
				const iterator it_f = end();
				while (it != it_f) {
					if (!new_ht.attempt_insertion(*it,tmp_key)) {
						// NOTE: here should we check with other hash functions before
						// giving up and increasing size?
						__PDEBUG(std::cout << "Cuckoo hash table resize triggered during resize.\n");
						const size_t new_index = new_ht.m_sizes_index + 1;
						__PDEBUG(std::cout << "Next size: " << sizes[new_index] << '\n');
						new_ht.resize(new_index);
						it = begin();
					} else {
						++it;
					}
				}
				swap(new_ht);
			}
			void attempt_shrink() {
				if (m_sizes_index == min_size_index) {
					__PDEBUG(std::cout << "Cuckoo hash table is already minimum size, won't shrink.\n");
					return;
				}
				cuckoo_hash_set new_ht;
				new_ht.resize(m_sizes_index - 1);
				// Assign current mults to the new table.
				new_ht.m_mults_index = m_mults_index;
				key_type tmp_key;
				const iterator it_f = end();
				for (iterator it = begin(); it != it_f; ++it) {
					if (!new_ht.attempt_insertion(*it,tmp_key)) {
						__PDEBUG(std::cout << "Unable to rebuild cuckoo hash set during shrinking, giving up.\n");
						return;
					}
				}
				swap(new_ht);
			}
			bool rehash() {
				cuckoo_hash_set new_ht;
				new_ht.m_mults_index = m_mults_index;
				new_ht.resize(m_sizes_index);
				for (size_t i = 0; i < max_rehash_tries; ++i) {
					key_type tmp_key;
					bool success = true;
					new_ht.next_mults();
					const iterator it_f = end();
					for (iterator it = begin(); it != it_f; ++it) {
						if (!new_ht.attempt_insertion(*it,tmp_key)) {
							__PDEBUG(std::cout << "Cuckoo hash table insertion failure during rehash, "
								"will try next mult.\n");
							success = false;
							break;
						}
					}
					if (success) {
						__PDEBUG(std::cout << "Rehash successful after " << (i + 1) << " tries.\n");
						// This means that we were able to insert all keys. Swap and return true.
						swap(new_ht);
						return true;
					}
					__PDEBUG(std::cout << "Rehash stopped at " << new_ht.size() << " out of " << size() << ".\n");
					// Let's reset new_ht.
					new_ht.resize(m_sizes_index);
				}
				__PDEBUG(std::cout << "Mults exhausted, rehash failed.\n")
				return false;
			}
			bool attempt_insertion(const key_type &k, key_type &tmp_key) {
				size_t pos = position1(hasher()(k));
				bucket_type *container = m_container;
				for (size_t i = 0; i < bucket_size; ++i) {
					// There's space in the bucket, rejoice!
					if (!container[pos].f[i]) {
						container[pos].f[i] = true;
						// NOTE: would it make sense here _not_ to default-initialise elements of the
						// bucket during construction and then use allocator::construct here instead of
						// assignment? Mhmh...
						container[pos].t[i] = k;
						++m_length;
						return true;
					}
				}
				// No space was found in the first-choice bucket. Choose randomly(?) the index of the element
				// in the bucket that will be displaced.
				const size_t dindex = pos & (bucket_size - 1);
				tmp_key = container[pos].t[dindex];
				container[pos].t[dindex] = k;
				size_t counter = 0;
				while (swap_and_displace(tmp_key,pos)) {
					++counter;
					if (counter > cuckoo_loop_threshold) {
						__PDEBUG(std::cout << "Cuckoo loop detected, returning false. Load factor is: " <<
							((static_cast<double>(m_length) + 1) / (sizes[m_sizes_index] * bucket_size)) << '\n');
						return false;
					}
				}
				return true;
			}
			// Place tmp_key into its location other than orig_location, displacing, if necessary, an existing
			// key. If displacement takes place, retval will be true and the content of tmp_key
			// will be the displaced one. Otherwise return false.
			bool swap_and_displace(key_type &tmp_key, size_t &orig_location) {
				// First thing we need to know if the original location was given by hash1 or hash2.
				const size_t h = hasher()(tmp_key), pos1 = position1(h),
					new_pos = (orig_location == pos1) ? position2(h) : pos1;
				bucket_type *container = m_container;
				for (size_t i = 0; i < bucket_size; ++i) {
					if (!container[new_pos].f[i]) {
						// Place found, rejoice!
						swapper()(container[new_pos].t[i],tmp_key);
						container[new_pos].f[i] = true;
						++m_length;
						return false;
					}
				}
				// No space was found in the alternative bucket. Choose randomly(?) the index of the element
				// in the bucket that will be displaced.
				const size_t dindex = new_pos & (bucket_size - 1);
				// Displace the selected key.
				container[new_pos].t[dindex].swap(tmp_key);
				orig_location = new_pos;
				return true;
			}
			static size_t m_hash(const size_t &h, const size_t &mult, const size_t &sizes_index) {
				return (mult * h) >> ((sizeof(max_fast_int) << 3) - sizes_index);
			}
			size_t position1(const size_t &h) const {
				p_assert(m_hash(h,mults[m_mults_index],m_sizes_index) < sizes[m_sizes_index]);
				return m_hash(h,mults[m_mults_index],m_sizes_index);
			}
			size_t position2(const size_t &h) const {
				p_assert(m_hash(h,mults[m_mults_index + 1],m_sizes_index) < sizes[m_sizes_index]);
				return m_hash(h,mults[m_mults_index + 1],m_sizes_index);
			}
			void next_mults() {
				m_mults_index += 2;
				if (m_mults_index == mults_size) {
					m_mults_index = 0;
				}
			}
		private:
			uint8				m_sizes_index;
			uint8				m_mults_index;
			size_t				m_length;
			bucket_type			*m_container;
	};

	template <class T, class HashFunction, class EqualKey, class Allocator, class SwapFunction>
	const size_t cuckoo_hash_set<T,HashFunction,EqualKey,Allocator,SwapFunction>::sizes[] = {
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

	// These are the decimal parts of the square roots of prime numbers multiplied by the
	// maximum value representable by size_t.
	// NOTE: here it is probably better to have all the double values below not below .5,
	// so that we are going to "randomise" almost all the bits of the result. Now these
	// work fine, so that it is probably better to wait for a stats framework of some kind
	// before attempting changes.
	template <class T, class HashFunction, class EqualKey, class Allocator, class SwapFunction>
	const size_t cuckoo_hash_set<T,HashFunction,EqualKey,Allocator,SwapFunction>::mults[] = {
		static_cast<size_t>(.7320508075688772 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.2360679774997898 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.6457513110645907 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.3166247903553998 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.6055512754639891 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.1231056256176606 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.3588989435406740 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.7958315233127191 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.3851648071345037 * boost::integer_traits<size_t>::const_max),
		static_cast<size_t>(.5677643628300215 * boost::integer_traits<size_t>::const_max)
	};
}

#endif
