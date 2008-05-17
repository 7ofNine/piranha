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

#include "base_series.h"
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
			void print(std::ostream &stream = std::cout, int limit = -1) const;
			std::string print_to_string() const;
			void save_to(const std::string &) const;
			void swap(Derived &);
			double norm() const;
			// TODO: maybe we can get rid of this with proper friends?
			const args_tuple_type &arguments() const;
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
			Derived partial(const std::string &) const;
			Derived partial(const psym &) const;
		protected:
			// TODO: check these protected methods, some of them can be moved into private
			// with proper friendship in manipulator classes.
			void construct_from_file(const std::string &);
			template <int N>
			void construct_from_psym(const psym &);
			void append_arg(const std::string &, const psym_p &);
			template <int N>
			void append_arg(const psym_p &);
			template <class Derived2>
			void merge_args(const Derived2 &);
			template <bool, class Derived2>
			Derived &merge_with_series(const Derived2 &);
			template <class Derived2>
			Derived &mult_by_series(const Derived2 &);
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

// Useful macros for named series.
#define NAMED_SERIES_TP_DECL class Cf, class Key, template <class> class I, \
				template <class, class, class, template <class> class> class Multiplier, \
				template <class> class Truncator, class Allocator
#define NAMED_SERIES_TP Cf,Key,I,Multiplier,Truncator,Allocator
#define REAL_NAMED_SERIES_TERM(term_name) term_name<Cf,Key,'|',Allocator>
#define REAL_NAMED_SERIES(series_name) series_name<NAMED_SERIES_TP>
#define REAL_NAMED_SERIES_BASE_ANCESTOR(term_name,series_name) piranha::base_series<REAL_NAMED_SERIES_TERM(term_name),'\n', \
	Allocator,REAL_NAMED_SERIES(series_name) >
#define REAL_NAMED_SERIES_NAMED_ANCESTOR(args,series_name) piranha::named_series<args,REAL_NAMED_SERIES(series_name) >

#define NAMED_SERIES_CTORS(series_name) \
	series_name() {nth_index<1>().max_load_factor(piranha::settings::load_factor());} \
	explicit series_name(const std::string &filename) \
	{ \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		named_ancestor::construct_from_file(filename); \
	} \
	explicit series_name(const piranha::max_fast_int &n) \
	{ \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		base_ancestor::construct_from_number(n,named_ancestor::m_arguments); \
	} \
	explicit series_name(const double &x) \
	{ \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		base_ancestor::construct_from_number(x,named_ancestor::m_arguments); \
	} \
	explicit series_name(const piranha::max_fast_int &n, const args_tuple_type &args_tuple) \
	{ \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
		base_ancestor::construct_from_number(n,args_tuple); \
	} \
	explicit series_name(const double &x, const args_tuple_type &args_tuple) \
	{ \
		nth_index<1>().max_load_factor(piranha::settings::load_factor()); \
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
