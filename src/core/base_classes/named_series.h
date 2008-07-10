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

#ifndef PIRANHA_NAMED_SERIES_H
#define PIRANHA_NAMED_SERIES_H

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../ntuple.h"
#include "../psym.h"
#include "../settings.h"

// Useful shortcuts.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_NAMED_SERIES_TP_DECL class ArgsDescr, class Derived
#define __PIRANHA_NAMED_SERIES_TP ArgsDescr,Derived

namespace piranha
{
	/// Named series toolbox.
	/**
	 * Toolbox for generating series with arguments.
	 * ArgsDescr must be a boost::tuple of structures each one containing a static const string
	 * called "name" naming the arguments of the series.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	class named_series
	{
			typedef ArgsDescr arguments_description;
			// Evaluation type. Used internally.
			typedef typename eval_type<Derived>::type eval_type;
		public:
			/// Compile-time constant for the number of arguments sets.
			static const int n_arguments_sets = boost::tuples::length<arguments_description>::value;
			BOOST_STATIC_ASSERT(n_arguments_sets > 0);
			typedef typename ntuple<vector_psym_p, n_arguments_sets>::type args_tuple_type;
			std::complex<Derived> complex() const;
			void print(std::ostream &stream = std::cout, int limit = -1) const;
			void save_to(const std::string &) const;
			template <class Cmp>
			void sort(const Cmp &);
			template <class Filter>
			Derived filter(const Filter &) const;
			void swap(Derived &);
			template <class Derived2>
			void merge_args(const Derived2 &);
			void trim();
			double norm() const;
			const args_tuple_type &arguments() const;
			void set_arguments(const args_tuple_type &);
			eval_type eval(const double &) const;
			Derived &operator+=(const max_fast_int &);
			Derived &operator+=(const double &);
			Derived &operator+=(const Derived &);
			Derived &operator-=(const max_fast_int &);
			Derived &operator-=(const double &);
			Derived &operator-=(const Derived &);
			Derived operator-() const;
			Derived &operator*=(const max_fast_int &);
			Derived &operator*=(const double &);
			Derived &operator*=(const Derived &);
			Derived &operator/=(const max_fast_int &);
			Derived &operator/=(const double &);
			Derived pow(const double &) const;
			Derived pow(const max_fast_int &) const;
			Derived root(const max_fast_int &) const;
			Derived partial(const std::string &) const;
			Derived partial(const psym &) const;
			template <class SubSeries>
			Derived sub(const psym &, const SubSeries &) const;
			template <class SubSeries>
			Derived sub(const std::string &, const SubSeries &) const;
			// Used in pyranha.
			std::string py_repr() const;
			template <class Term2>
			void py_append(const Term2 &);
			std::string py_arguments_description() const;
			args_tuple_type py_arguments() const;
			void py_shared_arguments_set() const;
		protected:
			// TODO: check these protected methods, some of them can be moved into private
			// with proper friendship in manipulator classes.
			void construct_from_file(const std::string &);
			template <int N>
			void construct_from_psym(const psym &);
			void append_arg(const std::string &, const psym_p &);
			template <int N>
			void append_arg(const psym_p &);
			template <bool, class Derived2>
			Derived &merge_with_series(const Derived2 &);
			template <class Derived2>
			Derived &mult_by_series(const Derived2 &);
			template <bool, class Number>
			Derived &merge_number_helper(const Number &);
			template <class Number>
			Derived &mult_number_helper(const Number &);
			template <class Number>
			Derived &divide_number_helper(const Number &);
		private:
			void print_plain(std::ostream &, int) const;
			void print_latex(std::ostream &, int) const;
			void read_from_file(std::ifstream &, const std::string &);
			void read_sections(std::ifstream &);
			void read_arg(std::ifstream &, const std::string &);
			void read_terms(std::ifstream &);
			template <class Derived2>
			bool is_args_compatible(const Derived2 &) const;
			template <class Derived2>
			void merge_incompatible_args(const Derived2 &);
			template <class Argument>
			Derived generic_partial(const Argument &) const;
			template <class Argument, class SubSeries>
			Derived generic_sub(const Argument &, const SubSeries &) const;
		protected:
			// Data members.
			args_tuple_type                 m_arguments;
			static std::vector<std::string> unknown_data;
	};

	// Initialization of static members.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	const int named_series<__PIRANHA_NAMED_SERIES_TP>::n_arguments_sets;

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	std::vector<std::string> named_series<__PIRANHA_NAMED_SERIES_TP>::unknown_data;

	// Meta-programming to get a tuple of (presence-flag + positional index) pairs for
	// a psym, given an arguments_tuple.
	template <class PosTuple, class ArgsTuple>
	struct named_series_get_psym_p_positions {
		p_static_check(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value, "");
		static void run(const psym_p &p, PosTuple &pos_tuple, const ArgsTuple &args_tuple) {
			// Set to not found.
			pos_tuple.template get_head().first = false;
			const size_t w = args_tuple.template get_head().size();
			for (size_t i = 0; i < w ; ++i) {
				if (args_tuple.template get_head()[i] == p) {
					pos_tuple.template get_head().first = true;
					pos_tuple.template get_head().second = i;
					break;
				}
			}
			named_series_get_psym_p_positions<typename PosTuple::tail_type, typename ArgsTuple::tail_type>::
			run(p, pos_tuple.template get_tail(), args_tuple.template get_tail());
		}
	};

	template <>
	struct named_series_get_psym_p_positions<boost::tuples::null_type, boost::tuples::null_type> {
		static void run(const psym_p &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

// Useful macros for named series.
#define E0_SERIES_NAMED_ANCESTOR(args,series_name) piranha::named_series<args,E0_SERIES(series_name) >

#define E1_SERIES_NAMED_ANCESTOR(args1,args2,series_name) piranha::named_series<boost::tuple<args1,args2>,series_name >

#define NAMED_SERIES_CTORS(series_name) \
	series_name() {} \
	explicit series_name(const std::string &filename) \
	{ \
		named_ancestor::construct_from_file(filename); \
	} \
	explicit series_name(const piranha::max_fast_int &n) \
	{ \
		base_ancestor::construct_from_number(n,named_ancestor::m_arguments); \
		named_ancestor::trim(); \
	} \
	explicit series_name(const double &x) \
	{ \
		base_ancestor::construct_from_number(x,named_ancestor::m_arguments); \
		named_ancestor::trim(); \
	} \
	explicit series_name(const piranha::max_fast_int &n, const args_tuple_type &args_tuple) \
	{ \
		base_ancestor::construct_from_number(n,args_tuple); \
	} \
	explicit series_name(const double &x, const args_tuple_type &args_tuple) \
	{ \
		base_ancestor::construct_from_number(x,args_tuple); \
	}
}

#include "named_series_io.h"
#include "named_series_manip.h"
#include "named_series_math.h"
#include "named_series_probe.h"

#undef derived_const_cast
#undef derived_cast
#undef __PIRANHA_NAMED_SERIES_TP_DECL
#undef __PIRANHA_NAMED_SERIES_TP

#endif
