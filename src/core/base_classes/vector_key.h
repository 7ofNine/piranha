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

#include <cstddef>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_same.hpp>

#include "../config.h"
#include "../exceptions.h"
#include "../memory.h"
#include "../Psym.h"


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
        static_assert(Position >= 0, "Wrong position.");

			// The index into the vector is the same as into the argDescr. This is how the the symbols get
		    // connected to the related numeric factors/exponents
			typedef std::vector<T, CountingAllocator<T, std::allocator<T> > > ContainerType;

		public:

			/// Type of contained data.
			typedef T value_type;
			/// Size type.
			typedef typename ContainerType::size_type size_type;
			/// Const iterator
			typedef typename ContainerType::const_iterator const_iterator;
			/// Iterator
			typedef typename ContainerType::iterator iterator;
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
			VectorKey(): container() {}

			/// Copy ctor.
			VectorKey(const VectorKey &other): container(other.container) {}


			/// Copy ctor, different position..
			template <int Position2, class Derived2>
			VectorKey(const VectorKey<T, Position2, Derived2> &other): container(other.container) {}


			/// Ctor from Psym.
			/**
			 * If the echelon level matches input integer n.
             * and the argstuple has one element and this element is the same as p given 
             * then push a 1 onto the container.
			 * 
             * this is so complicated. WHen is that construction actually used? Can we simplify it and make it more general?
             * Warning! there is no real verifyable relation from p to where the entry (1) is pushed onto the container!
			 */
			template <class ArgsTuple>
			VectorKey(const Psym &p, const int &n, const ArgsTuple &argsTuple): container()
			{
				(void)p;
				(void)argsTuple;
				// Construct only if the echelon level matches.
				if (n == Position) 
				{
					PIRANHA_ASSERT(argsTuple.template get<Position>().size() == 1 && argsTuple.template get<Position>()[0] == p);

					container.push_back(value_type(1));
				}
			}


			/// Swap content.
			void swap(VectorKey &other)
			{
				container.swap(other.container);
			}


			/// Is padding needed in order to be compatible with argsTuple?
			template <class ArgsTuple>
			bool needsPadding(const ArgsTuple &argsTuple) const
			{
				return (container.size() < argsTuple.template get<Position>().size());
			}


			/// Is this insertion-compatible with argsTuple?
			template <class ArgsTuple>
			bool isInsertable(const ArgsTuple &argsTuple) const
			{
				return (container.size() <= argsTuple.template get<Position>().size());
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
			void padRight(const ArgsTuple &argsTuple)
			{
				PIRANHA_ASSERT(argsTuple.template get<Position>().size() >= container.size());

				container.resize(boost::numeric_cast<size_type>(argsTuple.template get<Position>().size()));
			}


            //
			// Apply layout tuple.
            // change the layout of the key according to the layoutTuple element in Position
			// The position of the  elements in the container are reordered according to the corresponding layoutElement
			// A layout tuple is a tuple of vectors of pairs bool, std::size_t.
			//
			// ArgsTuple is not really used anywhere so why is it here.
			// the index into the layout determines where it goes in the modified key
			// the value in the pair of the layout (second) determines where it comes from (i.e. the source)
			// the boolean in the pair says if it is actually done at this index. If they are all false the new container contents will be an empty container
			// all values will be 0(or wahtever the initialization value of the used type is)
			//
			template <class LayoutTuple, class ArgsTuple>
			void applyLayout(LayoutTuple const &layoutTuple, ArgsTuple const &)
			{
                static_assert((boost::is_same<std::vector<std::pair<bool, std::size_t> >, typename boost::tuples::element<Position, LayoutTuple>::type>::value), "Wrong layout type.");
				// TODO: add check about tuples length.
				const size_type layoutSize = boost::numeric_cast<size_type>(layoutTuple.template get<Position>().size());

				// The layout must have at least all arguments in this.
				PIRANHA_ASSERT(layoutSize >= container.size());
				
                ContainerType newContainer(layoutSize);

				for (size_type i = 0; i < layoutSize; ++i) 
				{
					if (layoutTuple.template get<Position>()[i].first) // layout effective on *this
					{
						PIRANHA_ASSERT(layoutTuple.template get<Position>()[i].second < container.size());
						newContainer[i] = container[boost::numeric_cast<size_type>(layoutTuple.template get<Position>()[i].second)];
					}
				}

				newContainer.swap(container); // replace 
			}


            //
			// Test if vector key can be trimmed.
            // 
			template <class TrimFlags>
			void trimTest(TrimFlags &trimFlags) const
			{
				// TODO: add checks on TrimFlags type.
				const size_type size = container.size();
				
                PIRANHA_ASSERT(trimFlags.template get<Position>().size() == size);

				for (size_type i = 0; i < size; ++i) 
				{
					// If the element is different from zero, turn on the flag, don't trim it!
					if (container[i] != 0) 
					{
						trimFlags.template get<Position>()[i] = true;
					}
				}
			}


			/// Return trimmed version of this.
			// ArgsTuple is not used. Why is it here?
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &tf, const ArgsTuple &) const
			{
				// TODO: add checks on TrimFlags type.
				Derived retval;
				const size_type size = container.size();

				PIRANHA_ASSERT(tf.template get<position>().size() == size);

				// Make space, so we can avoid extra allocations in the cycle.
				retval.container.reserve(size);
				for (size_type i = 0; i < size; ++i) 
				{
					if (tf.template get<Position>()[i]) 
					{
						retval.container.push_back(container[i]);
					}
				}
				return retval;
			}


			/// Invert the sign of the integers in the array.
			void invertSign()
			{
				// NOTE: here perf with MP type could be improved by using in-place negate.
				const size_type size = container.size();
				for (size_type i = 0; i < size; ++i) 
				{
					container[i] = - container[i];
				}
			}


			/** @name Vector-like interface. */
			//@{
			/// Array-like operator[], const version.
			const value_type &operator[](const size_type &n) const
			{
				PIRANHA_ASSERT(n < container.size());
				return container[n];
			}


			/// Array-like operator[], mutable version.
			value_type &operator[](const size_type &n)
			{
				PIRANHA_ASSERT(n < container.size());

				return container[n];
			}


			/// Resize.
			void resize(const size_type &newSize)
			{
				container.resize(newSize);
			}


			/// Size.
			size_type size() const
			{
				return container.size();
			}


			/// Const begin.
			const_iterator begin() const
			{
				return container.begin();
			}


			/// Const end.
			const_iterator end() const
			{
				return container.end();
			}


			/// Begin.
			iterator begin()
			{
				return container.begin();
			}


			/// End.
			iterator end()
			{
				return container.end();
			}


			//@}
			/// Reverse lexicographic comparison.
			bool revlexComparison(const VectorKey &v2) const
			{
				const size_type size = this->size();
				PIRANHA_ASSERT(size == v2.size());

				// Shortcut in case there are no elements to compare.
				if (!size)
				{
					return false;
				}

				// Now we are certain that the size is at least 1, extract pointer to first element.
				// C++ standard guarantees that elements in std::vector are in contiguous memory areas.
				const value_type *ptr1 = &container[0];
                const value_type *ptr2 = &(v2.container[0]);
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
			bool lexComparison(const VectorKey &v2) const
			{
				const size_type size = this->size();
				PIRANHA_ASSERT(size == v2.size());
				if (!size) 
				{
					return false;
				}

				const value_type *ptr1 = &container[0];
                const value_type  *ptr2 = &(v2.container[0]);
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
				return (container == v2.container);
			}


			/// Equality test for elements.
			bool elementsEqualTo(const VectorKey &v2) const
			{
				return (container == v2.container);
			}


		protected:

			/// Print to stream the elements separated by the separator character.
			void printElements(std::ostream &outStream) const
			{
				const size_type size = this->size();
				for (size_type i = 0; i < size; ++i) 
				{
					outStream << container[i];
					// Print the separator iff this is not the last element.
					if (i != (size - 1))
                    {
						outStream << separator;
					}
				}
			}

			// the bool means we are actually printing it, the index in the pair mean which elemnt of the key we are actually printing
			// no separator is printed, the width is fixed to 6. we can improve here
			void printElementsSorted(std::ostream &outStream, std::vector<std::pair<bool, std::size_t> > positions) const
			{
				const size_type size = this->size();
				PIRANHA_ASSERT(size == positions.size());
				for (size_type i = 0; i < size; ++i)
				{
					if (positions[i].first)
					{
						outStream << std::setw(6) << std::right << container[positions[i].second] << " ";
					}
				}

			}

			/// Test for zero elements.
			/**
			 * Returns true if all elements are zero or size is zero, false otherwise.
			 */
			bool elementsAreZero() const
			{
				const size_type size = this->size();
				if (!size) 
				{
					return true;
				}

				const value_type *ptr = &container[0];
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
			std::size_t elementsHasher() const
			{
				const size_type size = this->size();

				if (!size) 
				{
					return 0;
				}

				std::size_t retval = 0;
				const value_type *ptr = &container[0];

				for (size_type i = 0; i < size; ++i) 
				{
					boost::hash_combine(retval, ptr[i]);
				}
				return retval;
			}

		protected:

			ContainerType container;
	};
};


#endif
