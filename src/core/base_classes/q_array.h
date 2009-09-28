/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redis\bute it and/or modify  *
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

#ifndef PIRANHA_Q_ARRAY_H
#define PIRANHA_Q_ARRAY_H

#include <algorithm>
#include <cstddef>
#include <functional>
#include <boost/functional/hash.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

#include "../config.h"
#include "../exceptions.h"
#include "../memory.h"
#include "../mp.h"
#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Array of rationals.
	/**
	 * Stores internally an array of piranha::mp_rational.
	 */
	template <int Pos, class Allocator, class Derived>
	class q_array {
			p_static_check(Pos >= 0, "Invalid position for q_array.");
		public:
			/// Make friends with q_array with other template parameters (i.e., position and/or allocator).
			template <int, class, class>
			friend class q_array;
			/// Size type.
			typedef uint8_t size_type;
			/// STL-like alias for mp_rational
			typedef mp_rational value_type;
			/// Allocator type.
			typedef counting_allocator<value_type,Allocator> allocator_type;
			/// Alias for position in arguments tuple.
			static const int position = Pos;
			/// Separator character for I/O.
			static const char separator = ';';
			/// Default ctor.
			/**
			 * Constructs an empty array.
			 */
			q_array(): m_size(0),m_ptr(0) {}
			/// Copy ctor.
			q_array(const q_array &other): m_size(other.m_size)
			{
				setup_elements_from_other(other);
			}
			/// Copy ctor. Position can be different.
			template <int Pos2, class Derived2>
			explicit q_array(const q_array<Pos2,Allocator,Derived2> &other): m_size(other.m_size)
			{
				setup_elements_from_other(other);
			}
			/// Ctor from psym.
			/**
			 * If the position matches input integer n, then
			 * resize to one element and set it to one.
			 */
			template <class ArgsTuple>
			q_array(const psym &p, const int &n, const ArgsTuple &args_tuple): m_size(0),m_ptr(0)
			{
				(void)p;
				(void)args_tuple;
				// Construct only if the positions match.
				if (n == position) {
					piranha_assert(args_tuple.template get<position>().size() == 1 &&
							 args_tuple.template get<position>()[0] == p);
					resize(1);
					m_ptr[0] = 1;
				}
			}
			/// Destructor.
			~q_array()
			{
				destroy_elements();
				deallocate_memory();
			}
			/// Assignment operator.
			q_array &operator=(const q_array &other)
			{
				if (this != &other) {
					// If the sizes are equal we can just copy over the elements,
					// without additional allocations.
					if (m_size == other.m_size) {
						std::copy(other.m_ptr, other.m_ptr + m_size, m_ptr);
					} else {
						destroy_elements();
						deallocate_memory();
						setup_elements_from_other(other);
					}
					// Take care of the size data member.
					m_size = other.m_size;
				}
				return *this;
			}
			/// Size getter.
			size_type size() const
			{
				return m_size;
			}
			/// Element getter.
			const value_type &operator[](const size_type &n) const
			{
				return m_ptr[n];
			}
			/// Swap with other q_array.
			void swap(q_array &other)
			{
				std::swap(m_size,other.m_size);
				std::swap(m_ptr,other.m_ptr);
			}
			/// Do I need padding in order to be compatible with args_tuple?
			template <class ArgsTuple>
			bool needs_padding(const ArgsTuple &args_tuple) const
			{
				return (m_size < args_tuple.template get<position>().size());
			}
			/// Am I insertion-compatible with args_tuple?
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &args_tuple) const
			{
				return (m_size <= args_tuple.template get<position>().size());
			}
			/// Number of atoms. Returns 1.
			std::size_t atoms() const
			{
				return 1;
			}
			/// Pad right.
			template <class ArgsTuple>
			void pad_right(const ArgsTuple &args_tuple)
			{
				piranha_assert(args_tuple.template get<position>().size() >= m_size);
				resize(boost::numeric_cast<size_type>(args_tuple.template get<position>().size()));
			}
			/// Apply layout.
			template <class Layout, class ArgsTuple>
			void apply_layout(const Layout &l, const ArgsTuple &)
			{
				const size_type l_size = boost::numeric_cast<size_type>(l.template get<position>().size());
				// The layout must have at least all arguments in this.
				piranha_assert(l_size >= m_size);
				q_array tmp;
				tmp.resize(l_size);
				for (size_type i = 0; i < l_size; ++i) {
					if (l.template get<position>()[i].first) {
						piranha_assert(l.template get<position>()[i].second < m_size);
						tmp[i] = (*this)[l.template get<Pos>()[i].second];
					}
				}
				swap(tmp);
			}
			/// Test need for trimming.
			template <class TrimFlags>
			void trim_test(TrimFlags &tf) const
			{
				piranha_assert(tf.template get<position>().size() == m_size);
				for (size_type i = 0; i < m_size; ++i) {
					// If the element is different from zero, turn on the flag..
					if ((*this)[i] != 0) {
						tf.template get<position>()[i] = true;
					}
				}
			}
			/// Return trimmed copy.
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &tf, const ArgsTuple &) const {
				Derived retval;
				std::vector<value_type> tmp;
				// Make space, so we can avoid extra allocations in the cycle.
				tmp.reserve(m_size);
				piranha_assert(tf.template get<position>().size() == m_size);
				for (size_type i = 0; i < m_size; ++i) {
					if (tf.template get<position>()[i]) {
						tmp.push_back((*this)[i]);
					}
				}
				retval.assign_vector(tmp);
				return retval;
			}
			/// Invert the sign of the integers in the array.
			void invert_sign() {
				for (size_type i = 0; i < m_size; ++i) {
					(*this)[i].negate();
				}
			}
			/// Upload to a std::vector.
			/**
			 * Vector size must not be smaller than current size, otherwise an assertion failure will be raised.
			 */
			void upload_to_vector(std::vector<value_type> &v) const
			{
				piranha_assert(v.size() >= m_size);
				for (size_type i = 0; i < m_size; ++i) {
					v[i] = (*this)[i];
				}
			}
			/// Assign a std::vector of values.
			template <class T>
			void assign_vector(const std::vector<T> &v)
			{
				const size_type size = boost::numeric_cast<size_type>(v.size());
				resize(size);
				for (std::size_t i = 0; i < size; ++i) {
					(*this)[i] = value_type(v[i]);
				}
			}
			/// Element-wise reverse lexicographic comparison.
			bool revlex_comparison(const Derived &a2) const
			{
				piranha_assert(m_size == a2.m_size);
				for (std::size_t i = m_size; i > 0; --i) {
					if ((*this)[i - 1] < a2[i - 1]) {
						return true;
					} else if ((*this)[i - 1] > a2[i - 1]) {
						return false;
					}
				}
				return false;
			}
			/// Element-wise lexicographic comparison.
			bool lex_comparison(const Derived &a2) const
			{
				piranha_assert(m_size == a2.m_size);
				for (std::size_t i = 0; i < m_size; ++i) {
					if ((*this)[i] < a2[i]) {
						return true;
					} else if ((*this)[i] > a2[i]) {
						return false;
					}
				}
				return false;
			}
		protected:
			/// Element setter.
			value_type &operator[](const size_type &n)
			{
				return m_ptr[n];
			}
			/// Print to stream the elements separated by the default separator character.
			void print_elements(std::ostream &out_stream) const
			{
				for (size_type i = 0; i < m_size; ++i) {
					out_stream << (*this)[i];
					// Print the separator iff this is not the last element.
					if (i != m_size - 1) {
						out_stream << separator;
					}
				}
			}
			/// Test for zero elements.
			/**
			 * Returns true if all elements are zero, false otherwise.
			 */
			bool elements_are_zero() const
			{
				for (size_type i = 0; i < m_size; ++i) {
					if ((*this)[i] != 0) {
						return false;
					}
				}
				return true;
			}
			/// Resize to new_size.
			/**
			 * Content will be unchanged to the minimum of old and new size.
			 */
			void resize(const size_type &new_size)
			{
				if (new_size == m_size) {
					return;
				}
				value_type *new_elements = 0;
				if (new_size > 0) {
					allocator_type a;
					new_elements = a.allocate(new_size);
					const size_type min_size = std::min(new_size,m_size);
					// Copy over old elements.
					for (size_type i = 0; i < min_size; ++i) {
						a.construct(new_elements + i, m_ptr[i]);
					}
					// Default-construct the remaining elements.
					const value_type tmp(0);
					for (size_type i = min_size; i < new_size; ++i) {
						a.construct(new_elements + i, tmp);
					}
				}
				// Destroy current elements.
				destroy_elements();
				deallocate_memory();
				// Assign data members.
				m_ptr = new_elements;
				m_size = new_size;
			}
			/// Hash value.
			std::size_t elements_hasher() const
			{
				std::size_t retval = 0;
				boost::hash<value_type> tmp_hash;
				for (size_type i = 0; i < m_size; ++i) {
					boost::hash_combine(retval, tmp_hash((*this)[i]));
				}
				return retval;
			}
			/// Equality test.
			bool elements_equal_to(const q_array &q) const
			{
				if (m_size == q.m_size) {
					for (size_type i = 0; i < m_size; ++i) {
						if ((*this)[i] != q[i]) {
							return false;
						}
					}
				} else {
					return false;
				}
				return true;
			}
		private:
			// Setup elements from other array: allocate the space needed to store
			// the elements and construct them from other. If other's size is 0, m_others
			// will be a null pointer.
			template <class T>
			void setup_elements_from_other(const T &other)
			{
				if (other.m_size > 0) {
					allocator_type a;
					m_ptr = a.allocate(other.m_size);
					for (size_type i = 0; i < other.m_size; ++i) {
						a.construct(m_ptr + i, other.m_ptr[i]);
					}
				} else {
					m_ptr = 0;
				}
			}
			// Destroy all allocated elements.
			void destroy_elements()
			{
				allocator_type a;
				for (size_type i = 0; i < m_size; ++i) {
					a.destroy(m_ptr + i);
				}
			}
			// Deallocate memory allocated in m_others.
			void deallocate_memory()
			{
				if (m_size > 0) {
					allocator_type().deallocate(m_ptr, m_size);
				}
			}
		private:
			size_type	m_size;
			value_type	*m_ptr;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
