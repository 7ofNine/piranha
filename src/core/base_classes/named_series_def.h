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

#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "../settings.h"
#include "../type_traits.h"
#include "base_series_def.h"
#include "named_series_mp.h"

#define __PIRANHA_NAMED_SERIES_TP_DECL class ArgsDescr, class Term, class Derived
#define __PIRANHA_NAMED_SERIES_TP            ArgsDescr, Term, Derived

namespace piranha
{
	/// Dictionary for evaluation.
	typedef boost::unordered_map<std::string, double> eval_dict;

	/// Named series toolbox.
	/**
	 * Toolbox for generating series with arguments.
	 * ArgsDescr must be a boost::tuple of structures each one containing a static const string
	 * called "name" naming the arguments of the series.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	class NamedSeries
	{
			template <class T, class Enable>
			friend struct named_series_add_selector;

			template <class T, class Enable>
			friend struct named_series_subtract_selector;

			template <class T, class Enable>
			friend struct named_series_multiply_selector;

			template <class T, class Enable>
			friend struct named_series_equality_selector;

		public:
			typedef ArgsDescr arguments_description;
			typedef typename Ntuple<VectorPsym, boost::tuples::length<arguments_description>::value>::type ArgsTupleType;

		private:
			struct s_iterator_generator
			{
				typedef Derived result_type;
				Derived operator()(const Term &t) const
				{
					Derived retval;
					retval.argumentsTuple = m_series.argumentsTuple;
					retval.insert(t, retval.argumentsTuple);
					return retval;
				}

				s_iterator_generator(const Derived &series):m_series(series) {}
				const Derived &m_series;
			};

		public:

			typedef boost::transform_iterator<s_iterator_generator, typename series_container<Term>::type::const_iterator> s_iterator;
			s_iterator s_begin() const;
			s_iterator s_end() const;
			std::complex<Derived> complex() const;

			void print(std::ostream &stream = std::cout) const;
			void print_plain(std::ostream &) const;
			void print_tex(std::ostream &) const;

			void save_to(const std::string &) const;
			// Rework this.
// 			template <class Filter>
// 			Derived filter(const Filter &) const;
			void swap(Derived &);
			double norm() const;
			typename term_eval_type_determiner<Term>::type eval(const double &) const;
			typename term_eval_type_determiner<Term>::type eval(const eval_dict &) const;
			
            std::size_t psi(const int &start = 0, const int &step = 1) const;
			
            const ArgsTupleType &arguments() const;
			
            void set_arguments(const ArgsTupleType &);

			template <class T>
			bool operator==(const T &) const;

			template <class T>
			bool operator!=(const T &) const;

			template <class T>
			Derived &operator+=(const T &);

			template <class T>
			Derived &operator-=(const T &);

			Derived operator-() const;

			template <class T>
			Derived &operator*=(const T &);

			template <class T>
			Derived &operator/=(const T &);

			Derived pow(const double &) const;
			Derived pow(const mp_rational &) const;
			Derived root(const int &) const;
			Derived partial(const std::string &, const int &n = 1) const;

			template <class SubSeries>
			Derived sub(const std::string &, const SubSeries &) const;

			template <class Key>
			Derived series_from_key(const Key &) const;

			template <class Cf>
			Derived series_from_cf(const Cf &) const;

			std::vector<std::vector<Derived> > split(const int &n = 0) const;

			std::vector<Derived> flatten() const;
			
            ~NamedSeries();
			
            //protected:
			void trim();
			template <class Derived2>
			void merge_args(const Derived2 &);
			// TODO: check these protected methods, some of them can be moved into private
			// with proper friendship in manipulator classes.
			void construct_from_file(const std::string &);
			template <int N>
			void construct_from_psym(const Psym &);

		private:

			template <class T>
			bool series_comparison(const T &) const;

			void append_arg(const std::string &, const Psym &);

			template <int N>
			void append_arg(const Psym &);

			template <class Derived2>
			Derived & mult_by_series(const Derived2 &);

			template <class Number>
			Derived & mult_number_helper(const Number &);

			template <class Number>
			Derived & divide_number_helper(const Number &);

			void print_pretty(std::ostream &) const;
			void read_from_file(std::ifstream &, const std::string &);
			void read_sections(std::ifstream &);
			void read_arg(std::ifstream &, const std::string &);
			void read_terms(std::ifstream &);

			template <class Derived2>
			bool is_args_compatible(const Derived2 &) const;

			template <class Derived2>
			void merge_incompatible_args(const Derived2 &);

			template <bool, class Derived2>
			Derived & merge_with_series(const Derived2 &);

		protected:

			// Data members.
			ArgsTupleType                   argumentsTuple;
			static std::vector<std::string> unknown_data;
	};

	// Initialization of static member.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	std::vector<std::string> NamedSeries<__PIRANHA_NAMED_SERIES_TP>::unknown_data;

// Useful macros for named series.
#define E0_SERIES_NAMED_ANCESTOR(args, termName, seriesName) piranha::NamedSeries<args, termName, E0_SERIES(seriesName)>

#define E1_SERIES_NAMED_ANCESTOR(args1, args2, termName, seriesName) piranha::NamedSeries<boost::tuple<args1, args2>, termName, seriesName>

#define NAMED_SERIES_BOILERPLATE(seriesName, N) \
public: \
	explicit seriesName() {} \
	explicit seriesName(const piranha::Psym &p) { \
		this->template construct_from_psym<N>(p); \
	} \
	explicit seriesName(const std::string &filename) \
	{ \
		this->construct_from_file(filename); \
	} \
	explicit seriesName(const double &x) \
	{ \
		*this = base_series_from_number(x, this->argumentsTuple); \
		this->trim(); \
	} \
	explicit seriesName(const piranha::mp_rational &q) \
	{ \
		*this = base_series_from_number(q, this->argumentsTuple); \
		this->trim(); \
	} \
	explicit seriesName(const piranha::mp_integer &z) \
	{ \
		*this = base_series_from_number(z, this->argumentsTuple); \
		this->trim(); \
	}
}

#endif
