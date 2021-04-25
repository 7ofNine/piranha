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

#ifndef PIRANHA_BASE_TERM_H
#define PIRANHA_BASE_TERM_H

#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <utility> // For std::pair.
#include <vector>

#include <boost/algorithm/string.hpp>

#include "../exceptions.h"
#include "../settings.h"

#define __PIRANHA_BASE_TERM_TP_DECL class Cf, class Key, char Separator, class Allocator, class Derived
#define __PIRANHA_BASE_TERM_TP Cf, Key, Separator, Allocator, Derived

namespace piranha
{
	template <int N, class Term>
	class BaseTermGetHelper
	{
		public:

			typedef typename Term::CfType Type;

			static Type &run(Term &t) 
            {
				return t.cf;
			}

			static const Type &run(const Term &t) 
            {
				return t.cf;
			}
	};

	template <class Term>
	class BaseTermGetHelper<1, Term>
	{
		public:

			typedef typename Term::KeyType Type;

			static Type &run(Term &t) 
            {
				return t.key;
			}

			static const Type &run(const Term &t) 
            {
				return t.key;
			}
	};


	// TODO: do we need this currently nowhere used!!
	/// Hasher functor.
	/**
	* Useful in STL-like containers.
	*/
	template<typename T>
	struct Hasher
	{
		std::size_t operator()(const T& t) const noexcept
		{
			return t.key.hash_value();
		}
	};


	// Base term class.
	//
	// Simple composition of coefficient and key classes.
	//
	// Cf:   coefficients for series term e.g. double_cf, polynomial_cf
	// Key:  key i.e. e.g. ExpoVector<boost::int16_t, 0>, TrigVector<boost::int16_t, 1>,
	//       last template parameter is actually the echelon level.
	// Separator: print/read separator between coefficiemt and key e.g.:  '|'
	// Allocator: specific allocator e.g. for statistics or performance improvements. but typicall std::allocator<char>
	// Derived: CRTP pattern, typically the derived class e.g. FourierSeriesTerm<Cf, Trig, Separator, Allocator>, or Monomial0

	// Concepts for Cf:  Cf(std::string, ArgsTuple) //constructor
	//                   Cf(cf2,         ArgsTuple) // constructor from another coefficient of a different kind)
	//                   Cf.swap(Cf)                // swap method 
	//                   Cf.print_plain(std::ostream,  ArgsTuple) // print method
	//                   Cf.print_pretty(std::ostream, ArgsTuple);
	//                   Cf.print_tex(std::ostream,    ArgsTuple)
	//                   operator==(int)            // equality to integer 
	//
	// Concepts for Key: Key(std::string, ArgsTuple) // constructor
	//                   Key(key2,       ArgsTuple)  // constructor from another key of a different kind) 
	//                   Key.swap(Key)               // swap method
	//                   Key.print_plain(std::ostream,  ArgsTuple) // print method
	//                   Key.print_pretty(std::ostream, ArgsTuple);
	//                   Key.print_tex(std::ostream,    ArgsTuple)
	//                   Key.is_unity()              // method Key value represents one. e.g. cos (0)
	//                   Key.operator==(Key)         // equality. Base terms only test the key
	//	                 Key.hash_value();

	template <__PIRANHA_BASE_TERM_TP_DECL>
	class BaseTerm
	{
		public:

			// Meta-programming to get the type of the component.
			// Is N actually used for anything?
			template <int N>
			struct Component 
			{
				typedef typename BaseTermGetHelper<N, BaseTerm>::Type Type;
			};

			/// Alias for coefficient type.
			typedef Cf CfType;
			/// Alias for key type.
			typedef Key KeyType;
			/// Alias for allocator type and interface
			using AllocatorType = typename std::allocator_traits<Allocator>::template rebind_alloc<Derived>;
			using AllocatorInterface = std::allocator_traits<AllocatorType>;




			/// Empty ctor.
			/**
			 * Default-initializes coefficient and key.
			 */
			BaseTerm(): cf(), key() {}
			
			// Ctor from string.
			// str is of type  "cf|key" for separator = "|"
			//What is ArgsTuple??
			template <class ArgsTuple>
			BaseTerm(const std::string &str, const ArgsTuple &argsTuple): cf(), key() 
            {
				std::vector<std::string> vs;
				boost::split(vs, str, boost::is_any_of(std::string(1, separator)));
				// should have precisely two components,
				if (vs.size() != 2) 
				{
					PIRANHA_THROW(value_error, std::string("unable to build term from input '") + str + "'");
				} else 
				{
					boost::trim(vs[0]); // coefficient
					boost::trim(vs[1]); // key
					// Try to build only if the strings actually contain something.
					if (!vs[0].empty()) 
                    {
						cf = CfType(vs[0], argsTuple);
					}
					if (!vs[1].empty()) 
                    {
						key = KeyType(vs[1], argsTuple);
					}
				}
			}
			

			// Copy ctor.
			//
			//Construct from BaseTerm with different coefficient and key.
		    // this requires a constuctor for cf(Derived2::CfType, ArgsTuple)
			//                            for key(Derived2::KeyType, ArgsTuple) 
			template <class Derived2, class ArgsTuple>
			BaseTerm(const Derived2 &t, const ArgsTuple &argsTuple)
			        : cf(t.cf, argsTuple), key(t.key) {}
			

			/// Ctor from coefficient - key pair.
			BaseTerm(const CfType &cf, const KeyType &key): cf(cf), key(key) {}

			template <int N>
			typename BaseTermGetHelper<N, BaseTerm>::Type &get() 
			{
				static_assert(N == 0 || N == 1);
				return BaseTermGetHelper<N, BaseTerm>::run(*this);
			}

			template <int N>
			const typename BaseTermGetHelper<N, BaseTerm>::Type &get() const 
			{
				static_assert(N == 0 || N == 1);
				return BaseTermGetHelper<N, BaseTerm>::run(*this);
			}

			void swap(BaseTerm &other) 
			{
				cf.swap(other.cf);
				key.swap(other.key);
			}

			// I/O.
			/// Print in plain format.
			template <class ArgsTuple>
			void printPlain(std::ostream &outStream, const ArgsTuple &argsTuple) const 
			{
				cf.printPlain(outStream, argsTuple); //print coefficient
				outStream << separator;
				key.printPlain(outStream, argsTuple); // print key
			}
			
			/// Print in pretty format.
			template <class ArgsTuple>
			void printPretty(std::ostream &outStream, const ArgsTuple &argsTuple) const 
			{
				if (key.isUnity()) 
				{
					cf.printPretty(outStream, argsTuple);
				} else if (cf == 1) 
				{
					key.printPretty(outStream, argsTuple);
				} else if (cf == -1) 
				{
					outStream << '-';
					key.printPretty(outStream, argsTuple);
				} else 
				{
					cf.printPretty(outStream, argsTuple);
					outStream << '*';
					key.printPretty(outStream, argsTuple);
				}
			}

			/// Print in tex format.
			template <class ArgsTuple>
			void printTex(std::ostream &outStream, const ArgsTuple &argsTuple) const 
			{
				if (key.isUnity()) 
				{
					cf.printTex(outStream, argsTuple);

				} else if (cf == 1) 
				{
					key.printTex(outStream, argsTuple);

				} else if (cf == -1) 
				{
					outStream << '-';
					key.printTex(outStream, argsTuple);

				} else 
				{
					cf.printTex(outStream, argsTuple);
					key.printTex(outStream, argsTuple);
				}
			}

			/// Equality test.
			/**
			 * Equality is defined by the equality of the keys.
			 */
			bool operator==(const BaseTerm &t) const 
			{
				return (key == t.key);
			}

			/// Check if the term is canonical.
            // TODO: why is this on the BAseTerm level???
			/**
			 * Will always return true, re-implement in derived term if necessary.
			 */
			template <class ArgsTuple>
			bool isCanonical(const ArgsTuple &) const 
			{
				return true;
			}

			/// Canonicalise the term.
            // TODO: why is this on the BAseTerm level???
			/**
			 * Won't do anything, re-implement in derived term if necessary.
			 */
			template <class ArgsTuple>
			void canonicalise(const ArgsTuple &) {}
			

			// TODO: do we need this currently nowhere used!!
			/// Hasher functor.
			/**
			 * Useful in STL-like containers.
			 */
			struct Hasher 
			{
				std::size_t operator()(const BaseTerm &t) const noexcept
				{
					return t.key.hash_value();
				}
			};

			// Data members.
			/// Coefficient.
			mutable CfType		cf;
			/// Key.
			KeyType		        key;
			/// Rebound allocator for term type.
			static AllocatorType	allocator;    //TODO: where is that actually used. Basseries uses is to allocate terms but is it necessary to prpagate it down to this level
			                                      // the allocator itself is not used in the baseTerm anywhere !!!
			/// Separator between coefficient and key in I/O.
			static const char separator = Separator;
	};

	// Static members initializations.
	template <__PIRANHA_BASE_TERM_TP_DECL>
	typename BaseTerm<__PIRANHA_BASE_TERM_TP>::AllocatorType
	BaseTerm<__PIRANHA_BASE_TERM_TP>::allocator;

	template <__PIRANHA_BASE_TERM_TP_DECL>
	const char BaseTerm<__PIRANHA_BASE_TERM_TP>::separator;

	/// Overload of hash_value function for piranha::BaseTerm.
    /// Don't rename is used by boost and name is required.
	/**
	 * The key's hash_value() method is used to calculate the term's hash value.
	 */
	template <__PIRANHA_BASE_TERM_TP_DECL>
	inline std::size_t hash_value(const BaseTerm<__PIRANHA_BASE_TERM_TP> &t)
	{
		return t.key.hash_value();
	}

#define PIRANHA_TERM_CTORS(TermName) \
	explicit TermName(): ancestor() {} \
	template <class ArgsTuple> \
	explicit TermName(const std::string &str, const ArgsTuple &argsTuple): \
			ancestor(str, argsTuple) {} \
	explicit TermName(const CfType &c, const KeyType &t): ancestor(c, t) {} \
	template <class Cf2, class ArgsTuple> \
	explicit TermName(const TermName<Cf2, KeyType, Separator, Allocator> &term, const ArgsTuple &a): \
			ancestor(term, a) {} \
	template <class Cf2, class Key2> \
	explicit TermName(const TermName<Cf2, Key2, Separator, Allocator> &term): \
			ancestor(term) {}
}

#undef __PIRANHA_BASE_TERM_TP_DECL
#undef __PIRANHA_BASE_TERM_TP

#endif
