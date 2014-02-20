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
#include <memory>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../memory.h"
#include "../psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Vector key.
	/**
	 * Series key type that can be represented as a vector of values.
	 */
	// T: type of key elements e.g boost::int16_t
	// Position: echelone level, determines which key it is. Each level has it's own key
	// Derived:  Derived class, for CRTP, static polymorphism
	template <class T, int Position, class Derived>
	class VectorKey
	{
			p_static_check(Position >= 0, "Wrong position.");

			typedef std::vector<T, counting_allocator<T, std::allocator<T> > > container_type;

		public:

			/// Type of contained data.
			typedef T value_type;
			/// Size type.
			typedef typename container_type::size_type size_type;
			/// Const iterator
			typedef typename container_type::const_iterator const_iterator;
			/// Iterator
			typedef typename container_type::iterator iterator;
			/// Position in the series' arguments tuple.
			// a series has Echelon+1 tuples as keys
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
			VectorKey(): m_container() {}

			/// Copy ctor.
			VectorKey(const VectorKey &other): m_container(other.m_container) {}


			/// Copy ctor, different position..
			template <int Position2, class Derived2>
			VectorKey(const VectorKey<T, Position2, Derived2> &other): m_container(other.m_container) {}


			/// Ctor from psym.
			/**
			 * If the position matches input integer n, then resize to one element and set it to one.
			 */
			template <class ArgsTuple>
			VectorKey(const psym &p, const int &n, const ArgsTuple &argsTuple): m_container()
			{
				(void)p;
				(void)argsTuple;
				// Construct only if the positions match.
				if (n == Position) 
				{
					piranha_assert(argsTuple.template get<Position>().size() == 1 && argsTuple.template get<Position>()[0] == p);
					m_container.push_back(value_type(1));
				}
			}


			/// Swap content.
			void swap(VectorKey &other)
			{
				m_container.swap(other.m_container);
			}


			/// Is padding needed in order to be compatible with argsTuple?
			template <class ArgsTuple>
			bool needs_padding(const ArgsTuple &argsTuple) const
			{
				return (m_container.size() < argsTuple.template get<Position>().size());
			}


			/// Is this insertion-compatible with argsTuple?
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &argsTuple) const
			{
				return (m_container.size() <= argsTuple.template get<Position>().size());
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
			void pad_right(const ArgsTuple &argsTuple)
			{
				piranha_assert(argsTuple.template get<Position>().size() >= m_container.size());

				m_container.resize(boost::numeric_cast<size_type>(argsTuple.template get<Position>().size()));
			}


			/// Apply layout tuple.
			/**
			 * A layout tuple is a tuple of vectors of pairs bool,std::size_t.
			 */
			template <class Layout, class ArgsTuple>
			void apply_layout(const Layout &l, const ArgsTuple &)
			{
				p_static_check((boost::is_same<std::vector<std::pair<bool, std::size_t> >, typename boost::tuples::element<Position, Layout>::type>::value), "Wrong layout type.");
				// TODO: add check about tuples length.
				const size_type l_size = boost::numeric_cast<size_type>(l.template get<Position>().size());

				// The layout must have at least all arguments in this.
				piranha_assert(l_size >= m_container.size());
				
                container_type new_container(l_size);

				for (size_type i = 0; i < l_size; ++i) 
				{
					if (l.template get<Position>()[i].first) 
					{
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

				for (size_type i = 0; i < size; ++i) 
				{
					// If the element is different from zero, turn on the flag.
					if (m_container[i] != 0) 
					{
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
				for (size_type i = 0; i < size; ++i) 
				{
					if (tf.template get<Position>()[i]) 
					{
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
				for (size_type i = 0; i < size; ++i) 
				{
					m_container[i] = - m_container[i];
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
			bool revlex_comparison(const VectorKey &v2) const
			{
				const size_type size = this->size();
				piranha_assert(size == v2.size());
				// Shortcut in case there are no elements to compare.
				if (!size)
				{
					return false;
				}
				// Now we are certain that the size is at least 1, extract pointer to first element.
				// C++ standard guarantees that elements in std::vector are in contiguous memory areas.
				const value_type *ptr1 = &m_container[0], *ptr2 = &(v2.m_container[0]);
				for (size_type i = size; i > 0; --i) 
				{
					if (ptr1[i - 1] < ptr2[i - 1]) 
					{
						return true;
					} else if (ptr1[i - 1] > ptr2[i - 1]) 
					{
						return false;
					}
				}
				return false;
			}


			/// Lexicographic comparison.
			bool lex_comparison(const VectorKey &v2) const
			{
				const size_type size = this->size();
				piranha_assert(size == v2.size());
				if (!size) 
				{
					return false;
				}
				const value_type *ptr1 = &m_container[0], *ptr2 = &(v2.m_container[0]);
				for (size_type i = 0; i < size; ++i)
				{
					if (ptr1[i] < ptr2[i]) 
					{
						return true;
					} else if (ptr1[i] > ptr2[i]) 
					{
						return false;
					}
				}
				return false;
			}


			/// Equality operator.
			bool operator==(const VectorKey &v2) const
			{
				return (m_container == v2.m_container);
			}


			/// Equality test for elements.
			bool elements_equal_to(const VectorKey &v2) const
			{
				return (m_container == v2.m_container);
			}


		protected:

			/// Print to stream the elements separated by the separator character.
			void print_elements(std::ostream &outStream) const
			{
				const size_type size = this->size();
				for (size_type i = 0; i < size; ++i) 
				{
					outStream << m_container[i];
					// Print the separator iff this is not the last element.
					if (i != (size - 1)) {
						outStream << separator;
					}
				}
			}


			/// Test for zero elements.
			/**
			 * Returns true if all elements are zero or size is zero, false otherwise.
			 */
			bool elements_are_zero() const
			{
				const size_type size = this->size();
				if (!size) 
				{
					return true;
				}

				const value_type *ptr = &m_container[0];
				for (size_type i = 0; i < size; ++i) 
				{
					if (ptr[i] != 0) 
					{
						return false;
					}
				}
				return true;
			}


			/// Hash value.
			/**
			 * Will produce a combined hash of all the elements of the vector using boost::hash_combine.
			 * An empty vector will produce a hash value of zero.
			 */
			std::size_t elements_hasher() const
			{
				const size_type size = this->size();
				if (!size) 
				{
					return 0;
				}
				std::size_t retval = 0;
				const value_type *ptr = &m_container[0];
				for (size_type i = 0; i < size; ++i) 
				{
					boost::hash_combine(retval, ptr[i]);
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
