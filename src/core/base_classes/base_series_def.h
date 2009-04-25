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

#include <boost/unordered_set.hpp>
#include <functional> // For std::equal_to.
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../memory.h"
#include "../psym.h"
#include "../type_traits.h"
#include "toolbox.h"

// Useful shortcuts.
#define __PIRANHA_BASE_SERIES_TP_DECL class Term, char Separator, class Allocator, class Derived
#define __PIRANHA_BASE_SERIES_TP Term,Separator,Allocator,Derived

namespace piranha
{
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	struct base_series {};

	/// Base series class.
	/**
	 * Term must derive from piranha::base_term class.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	class toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >
	{
			// Alias for term type.
			typedef Term term_type_;
			// Alias for coefficient type.
			typedef typename term_type_::cf_type cf_type;
			// Alias for key type.
			typedef typename term_type_::key_type key_type;
			// Alias for allocator type.
			typedef counting_allocator<term_type_,Allocator> allocator_type;
			// Term container.
			typedef boost::unordered_set<term_type_,boost::hash<term_type_>,std::equal_to<term_type_>,allocator_type>
				container_type;
		public:
			typedef term_type_ term_type;
			typedef typename term_eval_type_determiner<Term>::type base_eval_type;
			typedef typename container_type::iterator iterator;
			typedef typename container_type::const_iterator const_iterator;
			iterator begin();
			const_iterator begin() const;
			iterator end();
			const_iterator end() const;
			iterator find_term(const term_type &);
			const_iterator find_term(const term_type &) const;
			template <bool, bool, class Term2, class ArgsTuple>
			void insert(const Term2 &, const ArgsTuple &);
			template <class Term2, class ArgsTuple>
			void insert(const Term2 &, const ArgsTuple &);
			template <class Iterator, class ArgsTuple>
			void insert_range(const Iterator &, const Iterator &, const ArgsTuple &);
			template <class ArgsTuple>
			void term_erase(const iterator &, const ArgsTuple &);
			void clear_terms();
			void rehash(const size_t &);
			size_t length() const;
			bool empty() const;
			bool is_single_cf() const;
			size_t atoms() const;
			bool operator==(const Derived &) const;
			bool operator==(const double &) const;
			bool operator!=(const Derived &) const;
			bool operator!=(const double &) const;
		protected:
			template <class ArgsTuple>
			double base_norm(const ArgsTuple &) const;
			template <class ArgsTuple>
			base_eval_type base_eval(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived &base_add(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_add(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_subtract(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_subtract(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_mult_by(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_mult_by(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &base_divide_by(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived base_inv(const ArgsTuple &) const;
			template <class ArgsTuple>
			static Derived base_factorial(const int &, const ArgsTuple &);
			template <class ArgsTuple>
			static Derived base_choose(const int &, const int &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived base_pow(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived base_root(const int &, const ArgsTuple &) const;
			template <class PosTuple, class ArgsTuple>
			Derived base_partial(int, const PosTuple &, const ArgsTuple &) const;
			template <class PosTuple, class ArgsTuple>
			Derived base_partial(const PosTuple &, const ArgsTuple &) const;
			// Standard substitution functor. Will call sub() on coefficients and keys.
			struct sub_functor {
				template <class RetSeries, class Element, class PosTuple, class SubCaches,
					class ArgsTuple>
				static RetSeries run(const Element &e, const PosTuple &pos_tuple,
					SubCaches &sub_caches, const ArgsTuple &args_tuple) {
					return e.template sub<RetSeries>(pos_tuple, sub_caches, args_tuple);
				}
			};
			static const char separator = Separator;
			// Check that the separators do not conflict.
			p_static_check(separator != term_type::separator, "");
			template <class Number, class ArgsTuple>
			void construct_from_number(const Number &, const ArgsTuple &);
			template <class ArgsTuple>
			void base_construct_from_psym(const psym &, const int &, const ArgsTuple &);
			template <class ArgsTuple>
			void print_terms_plain(std::ostream &, const ArgsTuple &, int limit) const;
			template <class ArgsTuple>
			void print_terms_latex(std::ostream &, const ArgsTuple &, int limit) const;
			template <class ArgsTuple>
			void print_terms_pretty(std::ostream &, const ArgsTuple &, int limit) const;
			void swap_terms(Derived &);
			template <class Number>
			bool generic_numerical_comparison(const Number &) const;
			template <class ArgsTuple, class Layout>
			void apply_layout_to_terms(const ArgsTuple &, const Layout &, Derived &) const;
			template <bool, class Derived2, class ArgsTuple>
			Derived &merge_terms(const Derived2 &, const ArgsTuple &);
			template <class TrimFlags>
			void trim_test_terms(TrimFlags &) const;
			template <class TrimFlags, class ArgsTuple>
			void trim_terms(const TrimFlags &, Derived &, const ArgsTuple &) const;
			template <class T, class ArgsTuple>
			void multiply_coefficients_by(const T &, const ArgsTuple &);
			template <class T, class ArgsTuple>
			void divide_coefficients_by(const T &, const ArgsTuple &);
			template <bool, class Number, class ArgsTuple>
			Derived &merge_with_number(const Number &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived natural_power(const size_t &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived negative_integer_power(const int &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived real_power(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived nth_root(const int &, const ArgsTuple &) const;
			template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries base_sub(const PosTuple &, SubCaches &, const ArgsTuple &) const;
		private:
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
			template <class PosTuple, class ArgsTuple>
			static void ll_partial(const Derived &, Derived &, const PosTuple &, const ArgsTuple &);
			template <class ArgsTuple>
			bool common_power_handler(const double &, Derived &retval, const ArgsTuple &) const;
			template <class ArgsTuple>
			bool common_root_handler(const int &, Derived &retval, const ArgsTuple &) const;
		private:
			container_type m_container;
	};

#define E0_SERIES_TP_DECL class Cf, class Key, class Multiplier, class Truncator, class Allocator
#define E0_SERIES_TP Cf,Key,Multiplier,Truncator,Allocator
#define E0_SERIES_TERM(term_name) term_name<Cf,Key,'|',Allocator>
#define E0_SERIES(series_name) series_name<E0_SERIES_TP>
#define E0_SERIES_BASE_ANCESTOR(term_name,series_name) piranha::toolbox<piranha::base_series<E0_SERIES_TERM(term_name),'\n', \
	Allocator,E0_SERIES(series_name) > >

#define E1_SERIES_TP_DECL class Cf, class Key0, class Key1, \
						class Mult0, class Mult1, class Trunc0, class Trunc1, class Allocator
#define E1_SERIES_TP Cf,Key0,Key1,Mult0,Mult1,Trunc0,Trunc1,Allocator
#define E1_SERIES_COEFFICIENT(cf_name) cf_name<Cf,Key0,Mult0,Trunc0,Allocator>
#define E1_SERIES(series_name) series_name<E1_SERIES_TP>
#define E1_SERIES_TERM(term_name,cf_name) term_name< cf_name, Key1, '|', Allocator >
#define E1_SERIES_BASE_ANCESTOR(term_name,cf_name,series_name) piranha::toolbox<piranha::base_series<term_name< \
	cf_name,Key1,'|',Allocator>, \
	'\n',Allocator,series_name > >
}

#endif