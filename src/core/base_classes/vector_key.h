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

#ifndef PIRANHA_VECTOR_KEY_H
#define PIRANHA_VECTOR_KEY_H

#include <boost/functional/hash.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_same.hpp>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Vector key.
	/**
	 * Series key type that can be represented as a vector of values.
	 */
	template <class T, int Position, class Derived>
	class vector_key
	{
			p_static_check(Position >= 0,"Wrong position.");
			typedef std::vector<T> container_type;
		public:
			/// Type of contained data.
			typedef T value_type;
			/// Size type.
			typedef typename std::vector<T>::size_type size_type;
			/// Const iterator
			typedef typename std::vector<T>::const_iterator const_iterator;
			/// Iterator
			typedef typename std::vector<T>::iterator iterator;
			/// Position in the series' arguments tuple.
			static const int position = Position;
			/// Separator for string representation.
			/**
			 * The separator character must not be used in the textual representation of value_type,
			 * otherwise constructor from string will be confused.
			 */
			static const char separator = ';';
			/// Default ctor.
			/**
			 * Constructs an empty vector key..
			 */
			vector_key(): m_container() {}
			/// Copy ctor.
			vector_key(const vector_key &other): m_container(other.m_container) {}
			/// Copy ctor, different position..
			template <int Position2, class Derived2>
			explicit vector_key(const vector_key<T,Position2,Derived2> &other): m_container(other.m_container) {}
			/// Ctor from psym.
			/**
			 * If the position matches input integer n, then resize to one element and set it to one.
			 */
			template <class ArgsTuple>
			vector_key(const psym &p, const int &n, const ArgsTuple &args_tuple): m_container()
			{
				(void)p;
				(void)args_tuple;
				// Construct only if the positions match.
				if (n == Position) {
					piranha_assert(args_tuple.template get<Position>().size() == 1 && args_tuple.template get<Position>()[0] == p);
					m_container.push_back(value_type(1));
				}
			}
			/// Swap content.
			void swap(vector_key &other)
			{
				m_container.swap(other.m_container);
			}
			/// Is padding needed in order to be compatible with args_tuple?
			template <class ArgsTuple>
			bool needs_padding(const ArgsTuple &args_tuple) const
			{
				return (m_container.size() < args_tuple.template get<Position>().size());
			}
			/// Is this insertion-compatible with args_tuple?
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &args_tuple) const
			{
				return (m_container.size() <= args_tuple.template get<Position>().size());
			}
			/// Number of atoms.
			/**
			 * Will return 1.
			 */
			std::size_t atoms() const
			{
				return 1;
			}
			/// Pad right.
			template <class ArgsTuple>
			void pad_right(const ArgsTuple &args_tuple)
			{
				piranha_assert(args_tuple.template get<Position>().size() >= m_container.size());
				m_container.resize(boost::numeric_cast<size_type>(args_tuple.template get<Position>().size()));
			}
			/// Apply layout tuple.
			/**
			 * A layout tuple is a tuple of vectors of pairs bool,std::size_t.
			 */
			template <class Layout, class ArgsTuple>
			void apply_layout(const Layout &l, const ArgsTuple &)
			{
				p_static_check((boost::is_same<std::vector<std::pair<bool,std::size_t> >,typename boost::tuples::element<Position,Layout>::type>::value),"Wrong layout type.");
				// TODO: add check about tuples length.
				const size_type l_size = boost::numeric_cast<size_type>(l.template get<Position>().size());
				// The layout must have at least all arguments in this.
				piranha_assert(l_size >= m_container.size());
				container_type new_container(l_size);
				for (size_type i = 0; i < l_size; ++i) {
					if (l.template get<Position>()[i].first) {
						piranha_assert(l.template get<Position>()[i].second < m_container.size());
						new_container[i] = m_container[boost::numeric_cast<size_type>(l.template get<Position>()[i].second)];
					}
				}
				new_container.swap(m_container);
			}
			/// Test if vector key can be trimmed.
			template <class TrimFlags>
			void trim_test(TrimFlags &tf) const
			{
				// TODO: add checks on TrimFlags type.
				const size_type size = m_container.size();
				piranha_assert(tf.template get<Position>().size() == size);
				for (size_type i = 0; i < size; ++i) {
					// If the element is different from zero, turn on the flag.
					if (m_container[i] != 0) {
						tf.template get<Position>()[i] = true;
					}
				}
			}
			/// Return trimmed version of this.
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &tf, const ArgsTuple &) const
			{
				// TODO: add checks on TrimFlags type.
				Derived retval;
				const size_type size = m_container.size();
				piranha_assert(tf.template get<position>().size() == size);
				// Make space, so we can avoid extra allocations in the cycle.
				retval.m_container.reserve(size);
				for (size_type i = 0; i < size; ++i) {
					if (tf.template get<Position>()[i]) {
						retval.m_container.push_back(m_container[i]);
					}
				}
				return retval;
			}
			/// Invert the sign of the integers in the array.
			void invert_sign()
			{
				// NOTE: here perf with MP type could be improved by using in-place negate.
				const size_type size = m_container.size();
				for (size_type i = 0; i < size; ++i) {
					m_container[i] = -m_container[i];
				}
			}
			/** @name Vector-like interface. */
			//@{
			/// Array-like operator[], const version.
			const value_type &operator[](const size_type &n) const
			{
				piranha_assert(n < m_container.size());
				return m_container[n];
			}
			/// Array-like operator[], mutable version.
			value_type &operator[](const size_type &n)
			{
				piranha_assert(n < m_container.size());
				return m_container[n];
			}
			/// Resize.
			void resize(const size_type &new_size)
			{
				m_container.resize(new_size);
			}
			/// Size.
			size_type size() const
			{
				return m_container.size();
			}
			/// Const begin.
			const_iterator begin() const
			{
				return m_container.begin();
			}
			/// Const end.
			const_iterator end() const
			{
				return m_container.end();
			}
			/// Begin.
			iterator begin()
			{
				return m_container.begin();
			}
			/// End.
			iterator end()
			{
				return m_container.end();
			}
			//@}
			/// Reverse lexicographic comparison.
			bool revlex_comparison(const vector_key &v2) const
			{
				const size_type size = this->size();
				piranha_assert(size == v2.size());
				const value_type *ptr1 = &m_container[0], *ptr2 = &(v2.m_container[0]);
				for (size_type i = size; i > 0; --i) {
					if (ptr1[i - 1] < ptr2[i - 1]) {
						return true;
					} else if (ptr1[i - 1] > ptr2[i - 1]) {
						return false;
					}
				}
				return false;
			}
			/// Lexicographic comparison.
			bool lex_comparison(const vector_key &v2) const
			{
				const size_type size = this->size();
				piranha_assert(size == v2.size());
				const value_type *ptr1 = &m_container[0], *ptr2 = &(v2.m_container[0]);
				for (size_type i = 0; i < size; ++i) {
					if (ptr1[i] < ptr2[i]) {
						return true;
					} else if (ptr1[i] > ptr2[i]) {
						return false;
					}
				}
				return false;
			}
			/// Equality operator.
			bool operator==(const vector_key &v2) const
			{
				return (m_container == v2.m_container);
			}
			/// Equality test for elements.
			bool elements_equal_to(const vector_key &v2) const
			{
				return (m_container == v2.m_container);
			}
		protected:
			/// Print to stream the elements separated by the separator character.
			void print_elements(std::ostream &out_stream) const
			{
				const size_type size = this->size();
				for (size_type i = 0; i < size; ++i) {
					out_stream << m_container[i];
					// Print the separator iff this is not the last element.
					if (i != (size - 1)) {
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
				const value_type *ptr = &m_container[0];
				const size_type size = this->size();
				for (size_type i = 0; i < size; ++i) {
					if (ptr[i] != 0) {
						return false;
					}
				}
				return true;
			}
			/// Hash value.
			/**
			 * Will combine the hashes of all elements of the vector.
			 */
			std::size_t elements_hasher() const
			{
				const size_type size = this->size();
				const value_type *ptr = &m_container[0];
				std::size_t retval = 0;
				for (size_type i = 0; i < size; ++i) {
					boost::hash_combine(retval,ptr[i]);
				}
				return retval;
			}
		protected:
			container_type m_container;
	};
};

#undef derived_const_cast
#undef derived_cast

#endif
