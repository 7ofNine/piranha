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

#ifndef PIRANHA_BASE_SERIES_DEF_H
#define PIRANHA_BASE_SERIES_DEF_H

#include <boost/functional/hash.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/unordered_set.hpp>
#include <boost/utility/enable_if.hpp>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../memory.h"
#include "../mp.h"
#include "../psym.h"
#include "../type_traits.h"
#include "base_series_mp.h"
#include "base_series_tag.h"

// Template parameters list for piranha::base_series (declaration form).
#define __PIRANHA_BASE_SERIES_TP_DECL class Term, char Separator, class Allocator, class Derived
// Template parameters list for piranha::base_series (implementation form).
#define __PIRANHA_BASE_SERIES_TP Term,Separator,Allocator,Derived

namespace piranha
{
	// Implementation of echelon level determination.
	template <class CfSeries, class Enable = void>
	struct echelon_level_impl
	{
		static const int value = echelon_level_impl<typename CfSeries::term_type::cf_type>::value + 1;
	};

	template <class FinalCf>
	struct echelon_level_impl<FinalCf, typename boost::enable_if_c<!boost::is_base_of<base_series_tag,FinalCf>::value>::type>
	{
		static const int value = 0;
	};

	// Struct to define the container used in series - like a template typedef.
	template <class Term>
	struct series_container
	{
		typedef boost::unordered_set<Term,boost::hash<Term>,std::equal_to<Term>,counting_allocator<Term,std::allocator<char> > > type;
	};

	// These accessors are used in generic code that must work on both plain series (i.e., iterators) and sorted representations
	// of series as returned by get_sorted_series (i.e., pointers to pointers of terms).
	template <class Iterator>
	struct it_getter
	{
		static const Iterator &get(const Iterator &it)
		{
			return it;
		}
	};

	template <class TermPointer>
	struct it_getter<TermPointer *>
	{
		static const TermPointer get(const TermPointer *p)
		{
			return *p;
		}
	};

	/// Base series class.
	/**
	 * This class provides the basic representation of a series as a collection of terms stored into a hash set. The class is intended
	 * to be inherited together with (at least) either piranha::named_series (for a top-level series) or piranha::cf_series (for a coefficient series).
	 *
	 * The methods in this class provide the lowest level of series manipulation, and allow to operate directly on the individual terms of the
	 * series.
	 *
	 * @author Francesco Biscani (bluescarni@gmail.com)
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	class base_series: base_series_tag
	{
			// Befriend meta-programming classes.
			template <class RequestedCf, class SeriesCf>
			friend struct series_from_cf_impl;
			template <class RequestedKey, class SeriesKey>
			friend struct series_from_key_impl;
			template <int N>
			friend class series_flattener;
		public:
			/// Alias for term type.
			typedef Term term_type;
			/// Term container.
			/**
			 * The underlying term container is a plain boost::unordered set. Term types must specialise the boost::hash class, which
			 * will be used to provide the hash values for terms.
			 */
			typedef typename series_container<term_type>::type container_type;
			/// Echelon level.
			static const int echelon_level = echelon_level_impl<typename term_type::cf_type>::value;
			/// Type resulting from series evaluation.
			/**
			 * Determined automatically at compile time. Will be double in case of series with scalar coefficients,
			 * std::complex<double> in case of series with complex coefficients.
			 */
			typedef typename term_eval_type_determiner<Term>::type eval_type;
			/// Term separator in textual series representation.
			static const char separator = Separator;
			/// Const iterator over the series' terms.
			typedef typename container_type::const_iterator const_iterator;
			/// Size type.
			typedef typename container_type::size_type size_type;
			/** @name Series properties. */
			//@{
			size_type length() const;
			bool empty() const;
			bool is_single_cf() const;
			size_type atoms() const;
			//@}
		//protected:
			/** @name Construction/destruction. */
			//@{
			template <class Key, class ArgsTuple>
			static Derived base_series_from_key(const Key &, const ArgsTuple &);
			template <class Cf, class ArgsTuple>
			static Derived base_series_from_cf(const Cf &, const ArgsTuple &);
			template <class Number, class ArgsTuple>
			static Derived base_series_from_number(const Number &, const ArgsTuple &);
			template <class ArgsTuple>
			void base_construct_from_psym(const psym &, const int &, const ArgsTuple &);
			~base_series();
			//@}
			/** @name Series manipulation. */
			//@{
			void erase_term(const const_iterator &);
			void clear_terms();
			template <bool, bool, class Term2, class ArgsTuple>
			void insert(const Term2 &, const ArgsTuple &);
			template <class Term2, class ArgsTuple>
			void insert(const Term2 &, const ArgsTuple &);
			template <class Iterator, class ArgsTuple>
			void insert_range(const Iterator &, const Iterator &, const ArgsTuple &);
			void base_swap(Derived &);
			template <class Series, class ArgsTuple>
			void base_split(std::vector<std::vector<Series> > &, const int &n, const ArgsTuple &) const;
			template <class ArgsTuple>
			std::vector<term_type> flatten_terms(const ArgsTuple &) const;
			template <class Layout, class ArgsTuple>
			void apply_layout_to_terms(const Layout &, Derived &, const ArgsTuple &) const;
			template <bool, class Derived2, class ArgsTuple>
			Derived &merge_terms(const Derived2 &, const ArgsTuple &);
			template <class TrimFlags>
			void trim_test_terms(TrimFlags &) const;
			template <class TrimFlags, class ArgsTuple>
			void trim_terms(const TrimFlags &, Derived &, const ArgsTuple &) const;
			template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries base_sub(const PosTuple &, SubCaches &, const ArgsTuple &) const;
			template <bool, class Number, class ArgsTuple>
			Derived &merge_with_number(const Number &, const ArgsTuple &);
			//@}
			/** @name Terms accessors. */
			//@{
			const_iterator find_term(const term_type &) const;
			const_iterator begin() const;
			const_iterator end() const;
			//@}
			/** @name Base comparison methods. */
			//@{
			bool base_equal_to(const Derived &) const;
			bool base_equal_to(const double &) const;
			bool base_equal_to(const mp_rational &) const;
			bool base_equal_to(const mp_integer &) const;
			template <class Number>
			bool generic_numerical_comparison(const Number &) const;
			//@}
			/** @name Base maths. */
			//@{
			template <class T, class ArgsTuple>
			void multiply_coefficients_by(const T &, const ArgsTuple &);
			template <class T, class ArgsTuple>
			void divide_coefficients_by(const T &, const ArgsTuple &);
			template <class ArgsTuple>
			double base_norm(const ArgsTuple &) const;
			template <class ArgsTuple>
			eval_type base_eval(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived &base_add(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_add(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_add(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_add(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_subtract(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_subtract(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_subtract(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_subtract(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_mult_by(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_mult_by(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_mult_by(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_mult_by(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_divide_by(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_divide_by(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_divide_by(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived base_inv(const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived base_pow(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived base_pow(const mp_rational &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived natural_power(const std::size_t &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived negative_integer_power(const int &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived real_power(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived rational_power(const mp_rational &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived base_root(const int &, const ArgsTuple &) const;
			template <class Series, class PosTuple, class ArgsTuple>
			static void base_partial(const Derived &, Series &, const PosTuple &, const ArgsTuple &);
			template <class PosTuple, class ArgsTuple>
			Derived base_partial(int, const PosTuple &, const ArgsTuple &) const;
			template <class PosTuple, class ArgsTuple>
			Derived base_partial(const PosTuple &, const ArgsTuple &) const;
			//@}
			// Standard substitution functor. Will call sub() on coefficients and keys.
			struct sub_functor {
				template <class RetSeries, class Element, class PosTuple, class SubCaches, class ArgsTuple>
				static RetSeries run(const Element &e, const PosTuple &pos_tuple,
					SubCaches &sub_caches, const ArgsTuple &args_tuple)
				{
					return e.template sub<RetSeries>(pos_tuple, sub_caches, args_tuple);
				}
			};
			/** @name Base output streaming methods. */
			//@{
			template <class ArgsTuple>
			void print_terms_plain(std::ostream &, const ArgsTuple &) const;
			template <class ArgsTuple>
			void print_terms_tex(std::ostream &, const ArgsTuple &) const;
			template <class ArgsTuple>
			void print_terms_pretty(std::ostream &, const ArgsTuple &) const;
			//@}
		private:
			template <class Iterator, class ArgsTuple>
			void generic_print_terms_pretty(std::ostream &, const Iterator &, const Iterator &, const ArgsTuple &) const;
			template <class Iterator, class ArgsTuple>
			void generic_print_terms_tex(std::ostream &, const Iterator &, const Iterator &, const ArgsTuple &) const;
			template <class Iterator, class Series, class ArgsTuple>
			void generic_base_split(std::vector<std::vector<Series> > &, const Iterator &, const Iterator &, const ArgsTuple &) const;
			template <bool, class ArgsTuple>
			void ll_insert(const term_type &, const ArgsTuple &);
			template <bool, class ArgsTuple>
			void term_insert_new(const term_type &, const ArgsTuple &);
			template <int, class T, class ArgsTuple>
			void mult_div_coefficients_by(const T &, const ArgsTuple &);
			template <class Number, class ArgsTuple>
			Derived &multiply_by_number(const Number &, const ArgsTuple &);
			template <class Number, class ArgsTuple>
			Derived &divide_by_number(const Number &, const ArgsTuple &);
			template <class Number, class ArgsTuple>
			bool common_pow_handler(const Number &, Derived &retval, const ArgsTuple &) const;
		private:
			container_type m_container;
	};

#define E0_SERIES_TP_DECL class Cf, class Key, class Multiplier, class Truncator, class Allocator
#define E0_SERIES_TP Cf,Key,Multiplier,Truncator,Allocator
#define E0_SERIES_TERM(term_name) term_name<Cf,Key,'|',Allocator>
#define E0_SERIES(series_name) series_name<E0_SERIES_TP>
#define E0_SERIES_BASE_ANCESTOR(term_name,series_name) piranha::base_series<E0_SERIES_TERM(term_name),'\n', \
	Allocator,E0_SERIES(series_name) >

#define E1_SERIES_TP_DECL class Cf, class Key0, class Key1, \
						class Mult0, class Mult1, class Trunc0, class Trunc1, class Allocator
#define E1_SERIES_TP Cf,Key0,Key1,Mult0,Mult1,Trunc0,Trunc1,Allocator
#define E1_SERIES_COEFFICIENT(cf_name) cf_name<Cf,Key0,Mult0,Trunc0,Allocator>
#define E1_SERIES(series_name) series_name<E1_SERIES_TP>
#define E1_SERIES_TERM(term_name,cf_name) term_name< cf_name, Key1, '|', Allocator >
#define E1_SERIES_BASE_ANCESTOR(term_name,cf_name,series_name) piranha::base_series<term_name< \
	cf_name,Key1,'|',Allocator>, \
	'\n',Allocator,series_name >
}

#endif
