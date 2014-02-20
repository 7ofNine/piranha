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

#include <boost/algorithm/string.hpp>
#include <boost/static_assert.hpp>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <utility> // For std::pair.
#include <vector>

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

			typedef typename Term::cf_type type;

			static type &run(Term &t) 
            {
				return t.m_cf;
			}

			static const type &run(const Term &t) 
            {
				return t.m_cf;
			}
	};

	template <class Term>
	class BaseTermGetHelper<1, Term>
	{
		public:

			typedef typename Term::key_type type;

			static type &run(Term &t) 
            {
				return t.m_key;
			}

			static const type &run(const Term &t) 
            {
				return t.m_key;
			}
	};


	// Base term class.
	//
	// Simple composition of coefficient and key classes.
	//
	// Cf:   coefficients for series term e.g. double_cf, polynomial_cf
	// Key:  key i.e. e.g. TrigVector<boost::int16_t, 0>, TrigVector<boost::int16_t, 1>,
	//       last template parameter is actually the echelon level.
	// Separator: print/read separator between coefficiemt and key e.g.:  '|'
	// Allocator: specific allocator e.g. for statistics or performance improvements. but typicall std::allocator<char>
	// Derived: CRTP pattern, typically the derived class e.g. FourierSeriesTerm<Cf, Trig, Separator, Allocator>, 

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
			struct component 
			{
				typedef typename BaseTermGetHelper<N, BaseTerm>::type type;
			};

			/// Alias for coefficient type.
			typedef Cf cf_type;
			/// Alias for key type.
			typedef Key key_type;
			/// Alias for allocator type.
			typedef typename Allocator::template rebind<Derived>::other allocator_type;
			



			/// Empty ctor.
			/**
			 * Default-initializes coefficient and key.
			 */
			BaseTerm(): m_cf(), m_key() {}
			
			// Ctor from string.
			// str is of type  "cf|key" for separator = "|"
			//What is ArgsTuple??
			template <class ArgsTuple>
			BaseTerm(const std::string &str, const ArgsTuple &argsTuple): m_cf(), m_key() 
            {
				std::vector<std::string> vs;
				boost::split(vs, str, boost::is_any_of(std::string(1, separator)));
				// should have precisely two components,
				if (vs.size() != 2) 
				{
					piranha_throw(value_error, std::string("unable to build term from input '") + str + "'");
				} else 
				{
					boost::trim(vs[0]); // coefficient
					boost::trim(vs[1]); // key
					// Try to build only if the strings actually contain something.
					if (!vs[0].empty()) 
                    {
						m_cf = cf_type(vs[0], argsTuple);
					}
					if (!vs[1].empty()) 
                    {
						m_key = key_type(vs[1], argsTuple);
					}
				}
			}
			

			// Copy ctor.
			//
			//Construct from BaseTerm with different coefficient and key.
		    // this requires a constuctor for cf(Derived2::cf_type, ArgsTuple)
			//                            for key(Derived2::key_type, ArgsTuple) 
			template <class Derived2, class ArgsTuple>
			BaseTerm(const Derived2 &t, const ArgsTuple &argsTuple)
			        : m_cf(t.m_cf, argsTuple), m_key(t.m_key) {}
			

			/// Ctor from coefficient - key pair.
			BaseTerm(const cf_type &cf, const key_type &key): m_cf(cf), m_key(key) {}

			template <int N>
			typename BaseTermGetHelper<N, BaseTerm>::type &get() 
			{
				BOOST_STATIC_ASSERT(N == 0 or N == 1);
				return BaseTermGetHelper<N, BaseTerm>::run(*this);
			}

			template <int N>
			const typename BaseTermGetHelper<N, BaseTerm>::type &get() const 
			{
				BOOST_STATIC_ASSERT(N == 0 || N == 1);
				return BaseTermGetHelper<N, BaseTerm>::run(*this);
			}

			void swap(BaseTerm &other) 
			{
				m_cf.swap(other.m_cf);
				m_key.swap(other.m_key);
			}

			// I/O.
			/// Print in plain format.
			template <class ArgsTuple>
			void print_plain(std::ostream &outStream, const ArgsTuple &argsTuple) const 
			{
				m_cf.print_plain(outStream, argsTuple); //print coefficient
				outStream << separator;
				m_key.print_plain(outStream, argsTuple); // print key
			}
			
			/// Print in pretty format.
			template <class ArgsTuple>
			void print_pretty(std::ostream &outStream, const ArgsTuple &argsTuple) const 
			{
				if (m_key.is_unity()) 
				{
					m_cf.print_pretty(outStream, argsTuple);
				} else if (m_cf == 1) 
				{
					m_key.print_pretty(outStream, argsTuple);
				} else if (m_cf == -1) 
				{
					outStream << '-';
					m_key.print_pretty(outStream, argsTuple);
				} else 
				{
					m_cf.print_pretty(outStream, argsTuple);
					outStream << '*';
					m_key.print_pretty(outStream, argsTuple);
				}
			}

			/// Print in tex format.
			template <class ArgsTuple>
			void print_tex(std::ostream &outStream, const ArgsTuple &argsTuple) const 
			{
				if (m_key.is_unity()) 
				{
					m_cf.print_tex(outStream, argsTuple);

				} else if (m_cf == 1) 
				{
					m_key.print_tex(outStream, argsTuple);

				} else if (m_cf == -1) 
				{
					outStream << '-';
					m_key.print_tex(outStream, argsTuple);

				} else 
				{
					m_cf.print_tex(outStream, argsTuple);
					m_key.print_tex(outStream, argsTuple);
				}
			}

			/// Equality test.
			/**
			 * Equality is defined by the equality of the keys.
			 */
			bool operator==(const BaseTerm &t) const 
			{
				return (m_key == t.m_key);
			}

			/// Check if the term is canonical.
			/**
			 * Will always return true, re-implement in derived term if necessary.
			 */
			template <class ArgsTuple>
			bool is_canonical(const ArgsTuple &) const 
			{
				return true;
			}

			/// Canonicalise the term.
			/**
			 * Won't do anything, re-implement in derived term if necessary.
			 */
			template <class ArgsTuple>
			void canonicalise(const ArgsTuple &) {}
			
			/// Hasher functor.
			/**
			 * Useful in STL-like containers.
			 */
			struct hasher 
			{
				std::size_t operator()(const BaseTerm &t) const 
				{
					return t.m_key.hash_value();
				}
			};

			// Data members.
			/// Coefficient.
			mutable cf_type		m_cf;
			/// Key.
			key_type		    m_key;
			/// Rebound allocator for term type.
			static allocator_type	allocator;
			/// Separator between coefficient and key in I/O.
			static const char separator = Separator;
	};

	// Static members initializations.
	template <__PIRANHA_BASE_TERM_TP_DECL>
	typename BaseTerm<__PIRANHA_BASE_TERM_TP>::allocator_type
	BaseTerm<__PIRANHA_BASE_TERM_TP>::allocator;

	template <__PIRANHA_BASE_TERM_TP_DECL>
	const char BaseTerm<__PIRANHA_BASE_TERM_TP>::separator;

	/// Overload of hash_value function for piranha::BaseTerm.
	/**
	 * The key's hash_value() method is used to calculate the term's hash value.
	 */
	template <__PIRANHA_BASE_TERM_TP_DECL>
	inline std::size_t hash_value(const BaseTerm<__PIRANHA_BASE_TERM_TP> &t)
	{
		return t.m_key.hash_value();
	}

#define PIRANHA_TERM_CTORS(term_name) \
	explicit term_name(): ancestor() {} \
	template <class ArgsTuple> \
	explicit term_name(const std::string &str, const ArgsTuple &argsTuple): \
			ancestor(str, argsTuple) {} \
	explicit term_name(const cf_type &c, const key_type &t): ancestor(c, t) {} \
	template <class Cf2, class ArgsTuple> \
	explicit term_name(const term_name<Cf2, key_type, Separator, Allocator> &term, const ArgsTuple &a): \
			ancestor(term, a) {} \
	template <class Cf2, class Key2> \
	explicit term_name(const term_name<Cf2, Key2, Separator, Allocator> &term): \
			ancestor(term) {}
}

#undef __PIRANHA_BASE_TERM_TP_DECL
#undef __PIRANHA_BASE_TERM_TP

#endif
