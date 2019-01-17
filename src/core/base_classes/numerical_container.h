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

#ifndef PIRANHA_NUMERICAL_CONTAINER_H
#define PIRANHA_NUMERICAL_CONTAINER_H

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../math.h"
#include "../mp.h"
#include "../Psym.h"
#include "../settings.h"
#include "numerical_container_mp.h"
#include "numerical_container_tag.h"

// Convenience macros.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Numerical container class.
	/**
	 * This class can be used as a base class for coefficients that consist of a
	 * numerical entity (double, MP classes, etc.).
	 */
	template <class T, class Derived>
	class NumericalContainer: numerical_container_tag
	{
			template <class, class>
			friend struct numerical_container_constructor_selector;

			template <class, class, class>
			friend struct numerical_container_equality_selector;

			template <class, class>
			friend struct numerical_container_add_selector;

			template <class, class>
			friend struct numerical_container_subtract_selector;

			template <class, class>
			friend struct numerical_container_multiply_selector;

			template <class, class>
			friend struct numerical_container_divide_selector;

			template <class>
			friend struct numerical_container_print_tex_selector;

			template <class>
			friend struct numerical_container_print_pretty_selector;

		public:
            typedef Derived TermType;

			/// Typedef for evaluation type.
			typedef typename numerical_container_eval_type_determiner<T>::type EvalType;
			/// Alias for internal type.
			typedef T numerical_type;
			
			template <class, class SubCachesCons, class>
			struct SubstitutionCacheSelector
            {
				typedef SubCachesCons Type;
			};

			template <class, class SubCachesCons, class>
			struct EiSubstitutionCacheSelector
            {
				typedef SubCachesCons Type;
			};

			/// Default constructor, initialises internal value to 0.
			explicit NumericalContainer(): m_value(0) {}
			
			/// Constructor from string.
			/**
			 * Will call boost::lexical_converter internally.
			 */
			template <class ArgsTuple>
			explicit NumericalContainer(const std::string &s, const ArgsTuple &):
				m_value(boost::lexical_cast<T>(s)) {}
			
			/// Ctor from Psym.
			/**
			 * Sets internal value to one.
			 */
			template <class ArgsTuple>
			explicit NumericalContainer(const Psym &, const int &, const ArgsTuple &): m_value(1) {}
			
			/// Generic constructor.
			template <class U, class ArgsTuple>
			explicit NumericalContainer(const U &x, const ArgsTuple &):
				m_value(numerical_container_constructor_selector<U>::run(x))
			{}
			
			/// Print in plain mode.
			template <class ArgsTuple>
			void printPlain(std::ostream &outStream, const ArgsTuple &) const {
				outStream << boost::lexical_cast<std::string>(m_value);
			}
			
			/// Print in pretty mode. Equivalent to print_plain.
			template <class ArgsTuple>
			void printPretty(std::ostream &outStream, const ArgsTuple &) const {
				numerical_container_print_pretty_selector<T>::run(m_value,outStream);
			}
			
			/// Print in Tex mode.
			template <class ArgsTuple>
			void printTex(std::ostream &outStream, const ArgsTuple &) const 
			{
				numerical_container_print_tex_selector<T>::run(m_value, outStream);
			}
			
			/// Swap content using std::swap.
			Derived &swap(Derived &dc) 
			{
				std::swap(m_value, dc.m_value);
				return *derived_cast;
			}
			
			/// Pad right.
			template <class ArgsTuple>
			void padRight(const ArgsTuple &) {}
			
			/// Apply layout.
			template <class Layout, class ArgsTuple>
			void applyLayout(const Layout &, const ArgsTuple &) {}
			
			/// Test if trimming is possible.
			template <class TrimFlags>
			void trimTest(TrimFlags &) const {}
			
			/// Trim.
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &, const ArgsTuple &) const 
			{
				return *derived_const_cast;
			}
			
			/// Split.
			template <class Series, class ArgsTuple>
			void split(std::vector<std::vector<Series> > &, const int &, const ArgsTuple &) const
			{
				PIRANHA_ASSERT(false);
			}
			
			/// Number of atoms. Returns 1.
			std::size_t atoms() const {
				return 1;
			}

			/// Test for ignorability.
			/**
			 * Returns true if norm() is less than settings::get_numerical_zero().
			 */
			template <class ArgsTuple>
			bool isIgnorable(const ArgsTuple &a) const 
			{
				return (derived_const_cast->norm(a) < settings::get_numerical_zero());
			}
			
			/// Insertability test. Returns true.
			template <class ArgsTuple>
			bool isInsertable(const ArgsTuple &) const 
			{
				return true;
			}
			
			/// Padding test. Returns false.
			template <class ArgsTuple>
			bool needsPadding(const ArgsTuple &) const 
			{
				return false;
			}
			
			template <class U>
			bool operator==(const U &x) const 
			{
				return numerical_container_equality_selector<Derived,U>::run(*derived_const_cast,x);
			}

			template <class U>
			bool operator!=(const U &x) const 
			{
				return !operator==(x);
			}

			// Math.
			template <class ArgsTuple>
			void invertSign(const ArgsTuple &) 
			{
				m_value *= -1;
			}

			template <class U, class ArgsTuple>
			Derived &add(const U &x, const ArgsTuple &) 
			{
				return numerical_container_add_selector<U>::run(*derived_cast,x);
			}

			template <class U, class ArgsTuple>
			Derived &subtract(const U &x, const ArgsTuple &) 
			{
				return numerical_container_subtract_selector<U>::run(*derived_cast,x);
			}

			template <class U, class ArgsTuple>
			Derived &multBy(const U &x, const ArgsTuple &) 
			{
				return numerical_container_multiply_selector<U>::run(*derived_cast,x);
			}

			template <class U, class ArgsTuple>
			Derived &divideBy(const U &x, const ArgsTuple &) 
			{
				return numerical_container_divide_selector<U>::run(*derived_cast,x);
			}
			
			// Multiply and add.
			template <class Derived2, class ArgsTuple>
			void addmul(const Derived &x1, const Derived2 &x2, const ArgsTuple &) 
			{
				multiply_accumulate(m_value,x1.m_value,x2.get_value());
			}

			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &, const ArgsTuple &) const 
			{
				return Series();
			}
			
			template <class ArgsTuple>
			Derived pow(const double &x, const ArgsTuple &argsTuple) const 
			{
				return generic_pow(x,argsTuple);
			}

			template <class ArgsTuple>
			Derived pow(const mp_rational &q, const ArgsTuple &argsTuple) const 
			{
				return generic_pow(q,argsTuple);
			}
			
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &, SubCaches &, const ArgsTuple &argsTuple) const 
			{
				return RetSeries::baseSeriesFromCf(*derived_const_cast, argsTuple);
			}
			
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries eiSubstitute(const PosTuple &p, SubCaches &s, const ArgsTuple &a) const 
			{
				return sub<RetSeries>(p,s,a);
			}
			
			/// Bessel function of the first kind.
			/**
			 * Will use piranha::besselJ internally.
			 */
			template <class ArgsTuple>
			Derived besselJ(const int &n, const ArgsTuple &argsTuple) const
			{
				return Derived(piranha::besselJ(n,m_value),argsTuple);
			}
			
			/// Get value.
			const T &get_value() const 
			{
				return m_value;
			}
			
			/// Set value.
			template <class U>
			void set_value(const U &value) 
			{
				m_value = value;
			}

            /// get as string

            inline std::string toString() const
            {
                return boost::lexical_cast<std::string>(m_value);
            }

            inline std::size_t printLength() const
            {
                return toString().length();
            }

		private:

			template <class Number, class ArgsTuple>
			Derived generic_pow(const Number &n, const ArgsTuple &) const 
			{
				// COMPILER BUG: usual stuff with ubuntu's 4.1.3 compiler :(
				Derived retval;
				const numerical_type tmp(std::pow(m_value,n));
				retval.m_value = tmp;
				return retval;
			}

		private:

			T m_value;
	};


	// Overloads for I/O operators.
	template <class T, class Derived>
	inline std::istream &operator>>(std::istream &is, NumericalContainer<T, Derived> &nc)
	{
		std::string tmp;
		std::getline(is, tmp);
		nc = NumericalContainer<T, Derived>(boost::lexical_cast<T>(tmp));
		return is;
	}

	template <class T, class Derived>
	inline std::ostream &operator<<(std::ostream &os, const NumericalContainer<T, Derived> &nc)
	{
		os << boost::lexical_cast<std::string>(nc.get_value());
		return os;
	}

 //   template <class T, class Derived>
 //   inline std::size_t printLength(NumericalContainer<T, Derived> const &nc)
 //   {
 //       return boost::lexical_cast<std::string>(nc.get_value()).size();
 //   }

}

#undef derived_const_cast
#undef derived_cast

#endif
