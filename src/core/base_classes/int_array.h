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

#ifndef PIRANHA_INT_ARRAY_H
#define PIRANHA_INT_ARRAY_H

#include <algorithm> // For std::min/max.
#include <boost/integer.hpp>
#include <boost/integer_traits.hpp>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp> // For sub cache selection.
#include <stdint.h>
#include <utility> // For std::pair.
#include <vector>

#include "../config.h"
#include "../integer_typedefs.h"
#include "../math.h" // For lg.
#include "../memory.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Dynamically-sized integer array.
	/**
	 * Parametrized to a signed integer sized Bits, it contains also
	 * a flavour boolean flag, which can be used in trigonometric parts
	 * of Poisson series.
	 */
	template <int Bits, int Pos, class Allocator, class Derived>
	class int_array
	{
			p_static_check(Pos >= 0, "Invalid position for int_array.");
			p_static_check(Bits == 8 || Bits == 16, "Unsupported number of bits for value type in int_array.");
			typedef typename boost::int_t<Bits>::fast value_type_;
			typedef max_fast_int packed_type;
			p_static_check(sizeof(packed_type) % sizeof(value_type_) == 0,
				"Invalid packed/value ratio in int_array.");
			typedef counting_allocator<packed_type,Allocator> allocator_type;
			template <int Bits2, int Pos2, class Allocator2, class Derived2>
				friend class int_array;
			union container_type {
				value_type_ 	*v;
				packed_type		*p;
			};
			typedef uint8_t size_type_;
			static const size_type_ pack_capacity = sizeof(packed_type) / sizeof(value_type_);
			static const size_type_ pack_shift = lg<pack_capacity>::value;
		public:
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct sub_cache_selector {
				typedef boost::tuples::cons<typename Derived::template sub_cache<SubSeries,ArgsTuple>,
					SubCachesCons> type;
			};
			typedef value_type_ value_type;
			typedef size_type_ size_type;
			static const int position = Pos;
			static const char separator = ';';
			/// Default ctor.
			/**
			 * Constructs an empty array.
			 */
			int_array(): m_flavour(true), m_size(0) {
				m_container.p = allocator_type().allocate(0);
			}
			/// Copy ctor.
			int_array(const int_array &other): m_flavour(other.m_flavour), m_size(other.m_size) {
				init_copy(other);
			}
			/// Copy ctor. Position can be different.
			template <int Pos2, class Derived2>
			explicit int_array(const int_array<Bits,Pos2,Allocator,Derived2> &v): m_flavour(v.m_flavour),
				m_size(v.m_size) {
				init_copy(v);
			}
			/// Ctor from psym.
			/**
			 * If the position matches input integer n, then set the flavour to true,
			 * resize to one element and set it to one.
			 */
			template <class ArgsTuple>
			int_array(const psym_p &p, const int &n, const ArgsTuple &args_tuple):
					m_flavour(true), m_size(0) {
				(void)p;
				(void)args_tuple;
				m_container.p = allocator_type().allocate(0);
				// Construct only if the positions match.
				if (n == Pos) {
					p_assert(args_tuple.template get<Pos>().size() == 1 &&
							 args_tuple.template get<Pos>()[0] == p);
					resize(1);
					m_container.v[0] = 1;
				}
			}
			/// Dtor.
			~int_array() {
				allocator_type().deallocate(m_container.p, packed_size(m_size));
			}
			/// Assignment operator.
			int_array &operator=(const int_array &v) {
				// Don't do anything if the memory address is the same.
				if (this == &v) {
					return *this;
				}
				// Take care of flavour.
				m_flavour = v.m_flavour;
				const size_type p_size1 = packed_size(m_size), p_size2 = packed_size(v.m_size);
				if (p_size1 != p_size2) {
					allocator_type a;
					a.deallocate(m_container.p, p_size1);
					m_container.p = a.allocate(p_size2);
					m_size = v.m_size;
				}
				// Perform the copy from v to this.
				packed_copy(m_container.p, v.m_container.p, p_size2);
				return *this;
			}
			void swap(int_array &other) {
				std::swap(m_flavour,other.m_flavour);
				std::swap(m_size,other.m_size);
				std::swap(m_container.p,other.m_container.p);
			}
			/// Do I need padding in order to be compatible with args_tuple?
			template <class ArgsTuple>
			bool needs_padding(const ArgsTuple &args_tuple) const {
				return (m_size < args_tuple.template get<Pos>().size());
			}
			/// Am I insertion-compatible with args_tuple?
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &args_tuple) const {
				return (m_size <= args_tuple.template get<Pos>().size());
			}
			size_t atoms() const {
				return 1;
			}
			template <class ArgsTuple>
			bool checkup(const ArgsTuple &args_tuple) const {
				if (args_tuple.template get<Pos>().size() == m_size) {
					return true;
				}
				std::cout << "Size mismatch in int_array." << std::endl;
				return false;
			}
			/// Pad right.
			template <class ArgsTuple>
			void pad_right(const ArgsTuple &args_tuple) {
				p_assert(args_tuple.template get<Pos>().size() >= m_size);
				resize(args_tuple.template get<Pos>().size());
			}
			// TODO: invert order of parameters here.
			// TODO: check that l.get<Pos>().size() is compatible with size_type.
			template <class ArgsTuple, class Layout>
			void apply_layout(const ArgsTuple &, const Layout &l) {
				const size_t l_size = l.template get<Pos>().size();
				// The layout must have at least all arguments in this.
				p_assert(l_size >= m_size);
				container_type new_container;
				new_container.p = pack_init(l_size);
				for (size_t i = 0; i < l_size; ++i) {
					if (l.template get<Pos>()[i].first) {
						p_assert(l.template get<Pos>()[i].second < m_size);
						new_container.v[i] = m_container.v[l.template get<Pos>()[i].second];
					}
				}
				allocator_type().deallocate(m_container.p,packed_size(m_size));
				m_size = l_size;
				m_container.p = new_container.p;
			}
			template <class TrimFlags>
			void trim_test(TrimFlags &tf) const {
				p_assert(tf.template get<position>().size() == m_size);
				for (size_type i = 0; i < m_size; ++i) {
					// If the integer is different from zero, turn on the flag..
					if (m_container.v[i]) {
						tf.template get<position>()[i] = true;
					}
				}
			}
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &tf, const ArgsTuple &) const {
				Derived retval;
				retval.m_flavour = m_flavour;
				std::vector<int> tmp;
				p_assert(tf.template get<position>().size() == m_size);
				for (size_type i = 0; i < m_size; ++i) {
					if (tf.template get<position>()[i]) {
						tmp.push_back(m_container.v[i]);
					}
				}
				retval.assign_int_vector(tmp);
				return retval;
			}
			/// Invert the sign of the integers in the array.
			void invert_sign() {
				for (size_type i = 0; i < m_size; ++i) {
					m_container.v[i] = -m_container.v[i];
				}
			}
			template <class ArgsTuple>
			Derived inv_(const ArgsTuple &args_tuple) const {
				return derived_const_cast->pow_(-1,args_tuple);
			}
			/// Upload integers to vector of integers.
			void upload_ints_to(std::vector<int> &v) const {
				p_assert(v.size() >= m_size);
				for (size_type i = 0; i < m_size; ++i) {
					v[i] = m_container.v[i];
				}
			}
			/// Upload integers to vector of integer pairs.
			void upload_ints_to(std::vector<std::pair<int, int> > &v) const {
				p_assert(v.size() >= m_size);
				for (size_type i = 0; i < m_size; ++i) {
					v[i].first = m_container.v[i];
					v[i].second = m_container.v[i];
				}
			}
			/// Upload to v those integers which are less than the corresponding elements of v.
			void test_min_ints(std::vector<int> &v) const {
				p_assert(v.size() >= m_size);
				for (size_type i = 0; i < m_size; ++i) {
					if (m_container.v[i] < v[i]) {
						v[i] = m_container.v[i];
					}
				}
			}
			/// Upload to v.first/second those integers which are less than/greater than the corresponding elements of v.
			void test_min_max_ints(std::vector<std::pair<int, int> > &v) const {
				p_assert(v.size() >= m_size);
				for (size_type i = 0; i < m_size; ++i) {
					if (m_container.v[i] < v[i].first) {
						v[i].first = m_container.v[i];
					}
					if (m_container.v[i] > v[i].second) {
						v[i].second = m_container.v[i];
					}
				}
			}
			/// Codify integers into generalised lexicographic representation.
			template <class CodingVector, class ArgsTuple>
			max_fast_int code(const CodingVector &v, const ArgsTuple &) const {
				// The +1 is because the coding vector contains one extra element at the end.
				// The assert is >= instead of == beacuse we may code an array smaller than the
				// coding vector when multiplying series with different numbers of arguments.
				p_assert(v.size() >= static_cast<size_t>(m_size) + 1);
				max_fast_int retval = 0;
				for (size_type i = 0; i < m_size; ++i) {
					retval += (v[i] * m_container.v[i]);
				}
				return retval;
			}
			/// Decode integers from generalised lexicographic representation.
			template <class CodingVector, class MinMaxVec, class ArgsTuple>
			void decode(const max_fast_int &n, const CodingVector &cv, const max_fast_int &h_min,
						const MinMaxVec &mmv, const ArgsTuple &args_tuple) {
				resize(args_tuple.template get<position>().size());
				// The -1 is because the coding vector contains one extra element at the end.
				p_assert(cv.size() == static_cast<size_t>(m_size) + 1);
				const max_fast_int tmp = n - h_min;
				for (size_type i = 0; i < m_size; ++i) {
					m_container.v[i] = static_cast<value_type>((tmp % cv[i+1]) / cv[i] + mmv[i].first);
				}
			}
			void assign_int_vector(const std::vector<int> &v) {
				const size_t size = v.size();
				p_assert(boost::integer_traits<size_type>::const_max > size);
				// TODO: check where this function is used to see if this resize can be avoided.
				resize(size);
				// TODO: check for assignments out of numerical boundaries.
				for (size_t i = 0; i < size; ++i) {
					m_container.v[i] = (value_type)v[i];
				}
			}
			// Vector-like interface.
			/// Array-like operator[], const version.
			const value_type &operator[](const size_t &n) const {
				p_assert(n < m_size);
				return m_container.v[n];
			}
			/// Return container size.
			size_t size() const {
				return m_size;
			}
			bool revlex_comparison(const Derived &a2) const {
				p_assert(m_size == a2.m_size);
				value_type *ptr = m_container.v;
				for (size_t i = m_size; i > 0; --i) {
					if (ptr[i - 1] < a2[i - 1]) {
						return true;
					} else if (ptr[i - 1] > a2[i - 1]) {
						return false;
					}
				}
				return false;
			}
			bool lex_comparison(const Derived &a2) const {
				p_assert(m_size == a2.m_size);
				for (size_t i = 0; i < m_size; ++i) {
					if (m_container.v[i] < a2[i]) {
						return true;
					} else if (m_container.v[i] > a2[i]) {
						return false;
					}
				}
				return false;
			}
		protected:
			/// Array-like operator[], mutable version.
			value_type &operator[](const size_t &n) {
				p_assert(n < m_size);
				return m_container.v[n];
			}
			/// Print to stream the elements separated by the default separator character.
			void print_elements(std::ostream &out_stream) const {
				for (size_t i = 0; i < m_size; ++i) {
					// We cast to max_fast_int, which is the largest integer type admitted.
					out_stream << static_cast<int>(m_container.v[i]);
					// Print the separator iff this is not the last element.
					if (i != static_cast<size_t>(m_size - 1)) {
						out_stream << separator;
					}
				}
			}
			/// Test for zero elements.
			/**
			 * Returns true if all elements are zero, false otherwise.
			 */
			bool elements_are_zero() const {
				const size_type p_size = packed_size(m_size);
				for (size_type i = 0; i < p_size; ++i) {
					if (m_container.p[i] != 0) {
						return false;
					}
				}
				return true;
			}
			// TODO: should we pass size_t here, and test against size_type and throw in case of out-of-range request?
			/// Resize the container.
			/**
			 * The existing elements are copied over.
			 */
			void resize(const size_type &new_size) {
				if (m_size == new_size) {
					return;
				}
				container_type new_container;
				new_container.p = pack_init(new_size);
				// Copy to the minimum of the new sizes.
				const size_type old_p_size = packed_size(m_size);
				packed_copy(new_container.p, m_container.p,
					std::min<size_type>(old_p_size, packed_size(new_size)));
				// Destroy & assign new.
				allocator_type().deallocate(m_container.p, old_p_size);
				m_container.p = new_container.p;
				m_size = new_size;
			}
			/// Hash value.
			/**
			 * Hashes only the integer elements of the array, not the flavour.
			 */
			size_t elements_hasher() const {
				const size_type p_size = packed_size(m_size);
				switch (p_size) {
					case 0:
						return 0;
					case 1:
						return m_container.p[0];
				}
				size_t retval = m_container.p[0];
				for (size_type i = 1; i < p_size; ++i) {
					boost::hash_combine(retval, m_container.p[i]);
				}
				return retval;
			}
			/// Equality test.
			/**
			 * Tests only the integer elements of the array, not the flavour.
			 */
			bool elements_equal_to(const int_array &v) const {
				if (m_size == v.m_size) {
					const size_type p_size = packed_size(m_size);
					for (size_type i = 0; i < p_size; ++i) {
						if (m_container.p[i] != v.m_container.p[i]) {
							return false;
						}
					}
				} else {
					return false;
				}
				return true;
			}
		private:
			static size_type packed_size(const size_type &s) {
				return ((s >> pack_shift) + ((s & (pack_capacity - 1)) != 0));
			}
			static void packed_copy(packed_type *dest, const packed_type *src, const size_type &packed_size) {
				for (size_type i = 0; i < packed_size; ++i) {
					dest[i] = src[i];
				}
			}
			static packed_type *pack_init(const size_type &l_size) {
				const size_type p_size = packed_size(l_size);
				packed_type *retval = allocator_type().allocate(p_size);
				for (size_type i = 0; i < p_size; ++i) {
					retval[i] = 0;
				}
				return retval;
			}
			template <class IntArray>
			void init_copy(const IntArray &other) {
				const size_type p_size = packed_size(m_size);
				m_container.p = allocator_type().allocate(p_size);
				packed_copy(m_container.p, other.m_container.p, p_size);
			}
		protected:
			// Data members.
			/// Flavour.
			bool					m_flavour;
			/// Size of the array.
			/**
			 * Equal to the number of elements contained by the array.
			 */
			size_type				m_size;
			/// Container.
			container_type			m_container;
	};
};

#undef derived_const_cast
#undef derived_cast

#endif
