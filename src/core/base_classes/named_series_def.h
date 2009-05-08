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

#ifndef PIRANHA_NAMED_SERIES_DEF_H
#define PIRANHA_NAMED_SERIES_DEF_H

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../ntuple.h"
#include "../psym.h"
#include "../settings.h"
#include "../type_traits.h"
#include "toolbox.h"

#define __PIRANHA_NAMED_SERIES_TP_DECL class ArgsDescr, class Term, class Derived
#define __PIRANHA_NAMED_SERIES_TP ArgsDescr,Term,Derived

namespace piranha
{
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	struct named_series {};

	/// Named series toolbox.
	/**
	 * Toolbox for generating series with arguments.
	 * ArgsDescr must be a boost::tuple of structures each one containing a static const string
	 * called "name" naming the arguments of the series.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	class toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >
	{
		public:
			typedef ArgsDescr arguments_description;
			/// Compile-time constant for the number of arguments sets.
			static const int n_arguments_sets = boost::tuples::length<arguments_description>::value;
			p_static_check(n_arguments_sets > 0, "The number of arguments vector must be strictly positive.");
			typedef typename ntuple<vector_psym, n_arguments_sets>::type args_tuple_type;
			typedef typename term_eval_type_determiner<Term>::type eval_type;
			std::complex<Derived> complex() const;
			void print(std::ostream &stream = std::cout, int limit = -1) const;
			void save_to(const std::string &) const;
			// Rework this.
// 			template <class Filter>
// 			Derived filter(const Filter &) const;
			void swap(Derived &);
			double norm() const;
			eval_type eval(const double &) const;
			size_t psi(const int &start = 0, const int &step = 1) const;
			const args_tuple_type &arguments() const;
			void set_arguments(const args_tuple_type &);
			bool operator==(const Derived &) const;
			bool operator!=(const Derived &) const;
			bool operator==(const double &) const;
			bool operator!=(const double &) const;
			Derived &operator+=(const double &);
			Derived &operator+=(const Derived &);
			Derived &operator-=(const double &);
			Derived &operator-=(const Derived &);
			Derived operator-() const;
			Derived &operator*=(const double &);
			Derived &operator*=(const Derived &);
			Derived &operator/=(const double &);
			static Derived factorial(const int &);
			static Derived choose(const int &, const int &);
			Derived pow(const double &) const;
			Derived root(const int &) const;
			Derived inv() const;
			Derived partial(const psym &, const int &n = 1) const;
			template <class SubSeries>
			Derived sub(const psym &, const SubSeries &) const;
		protected:
			template <class Derived2>
			void merge_args(const Derived2 &);
			void trim();
			// TODO: check these protected methods, some of them can be moved into private
			// with proper friendship in manipulator classes.
			void construct_from_file(const std::string &);
			template <int N>
			void construct_from_psym(const psym &);
			void append_arg(const std::string &, const psym &);
			template <int N>
			void append_arg(const psym &);
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
			void print_pretty(std::ostream &, int) const;
			void print_latex(std::ostream &, int) const;
			void read_from_file(std::ifstream &, const std::string &);
			void read_sections(std::ifstream &);
			void read_arg(std::ifstream &, const std::string &);
			void read_terms(std::ifstream &);
			template <class Derived2>
			bool is_args_compatible(const Derived2 &) const;
			template <class Derived2>
			void merge_incompatible_args(const Derived2 &);
		protected:
			// Data members.
			args_tuple_type                 m_arguments;
			static std::vector<std::string> unknown_data;
	};

	// Initialization of static member.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	std::vector<std::string> toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::unknown_data;

// Useful macros for named series.
#define E0_SERIES_NAMED_ANCESTOR(args, term_name, series_name) piranha::toolbox<piranha::named_series<args, term_name, E0_SERIES(series_name) > >

#define E1_SERIES_NAMED_ANCESTOR(args1,args2, term_name, series_name) piranha::toolbox<piranha::named_series<boost::tuple<args1,args2>,term_name,series_name > >

#define NAMED_SERIES_BOILERPLATE(series_name,N) \
public: \
	explicit series_name() {} \
	explicit series_name(const piranha::psym &p) { \
		this->template construct_from_psym<N>(p); \
	} \
	explicit series_name(const std::string &filename) \
	{ \
		this->construct_from_file(filename); \
	} \
	explicit series_name(const double &x) \
	{ \
		this->construct_from_number(x,this->m_arguments); \
		this->trim(); \
	}
}

#endif