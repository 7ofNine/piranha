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

#ifndef PIRANHA_DEGREE_TRUNCATOR_H
#define PIRANHA_DEGREE_TRUNCATOR_H

#include <algorithm> // For sorting.
#include <boost/tuple/tuple.hpp>
#include <cmath> // For std::ceil.
#include <iostream>
#include <utility>
#include <vector>

#include "../base_classes/toolbox.h"
#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../psym.h"
#include "../settings.h" // For debug messages.

namespace piranha
{
	struct degree_truncator_ {};

	/// Truncator based on the minium degree of the series.
	template <>
	class __PIRANHA_VISIBLE toolbox<degree_truncator_>
	{
			enum mode {
				deg,
				p_deg,
				inactive
			};
		public:
			template <class Multiplier>
			class get_type
			{
					typedef typename ntuple<std::vector<std::pair<bool,size_t> >,
						boost::tuples::length<typename Multiplier::args_tuple_type>::value>::type
						pos_tuple_type;
					static const int expo_term_pos = Multiplier::series_type1::expo_term_position;
					static const int expo_args_pos = Multiplier::series_type1::expo_args_position;
					struct min_degree_comparison
					{
						template <class Term>
						bool operator()(const Term *t1, const Term *t2) const
						{
							return t1->template get<expo_term_pos>().min_degree() <
								   t2->template get<expo_term_pos>().min_degree();
						}
					};
					struct partial_min_degree_comparison
					{
						partial_min_degree_comparison(const pos_tuple_type &pos_tuple):m_pos_tuple(pos_tuple) {}
						template <class Term>
						bool operator()(const Term *t1, const Term *t2) const
						{
							return t1->template get<expo_term_pos>().partial_min_degree(m_pos_tuple) <
								   t2->template get<expo_term_pos>().partial_min_degree(m_pos_tuple);
						}
						const pos_tuple_type &m_pos_tuple;
					};
				public:
					typedef get_type type;
					get_type(Multiplier &m, bool initialise = true): m_multiplier(m)
					{
						// Some static checks.
						p_static_check(Multiplier::series_type1::expo_args_position ==
									   Multiplier::series_type2::expo_args_position, "");
						p_static_check(Multiplier::series_type1::expo_term_position ==
									   Multiplier::series_type2::expo_term_position, "");
						// Convert psyms vector into position tuple only if we are truncating to partial degree.
						if (m_mode == p_deg) {
							m_pos_tuple = psyms2pos(m_psyms,m_multiplier.m_args_tuple);
						}
						if (initialise) {
							init();
						}
					}
					template <class T>
					bool accept(const T &) const
					{
						return true;
					}
					template <class Term1, class Term2>
					bool skip(const Term1 &t1, const Term2 &t2) const
					{
						switch (m_mode) {
							case deg:
								return (t1.template get<expo_term_pos>().min_degree() +
									t2.template get<expo_term_pos>().min_degree() >=
									m_degree_limit);
							case p_deg:
								return (t1.template get<expo_term_pos>().partial_min_degree(m_pos_tuple) +
									t2.template get<expo_term_pos>().partial_min_degree(m_pos_tuple) >=
									m_degree_limit);
							case inactive:
								// We should never get there.
								piranha_assert(false);
						}
						return false;
					}
					// Number of a iterations of a power series development of a power series.
					// NOTE: if start is negative, it is assumed that negative powers of the input series
					// have a minimum degree which is proportional to the input series' and with its sign changed.
					template <class PowerSeries, class ArgsTuple>
					static size_t power_series_iterations(const PowerSeries &s, const int &start, const int &step_size,
						const ArgsTuple &args_tuple)
					{
						if (step_size < 1) {
							piranha_throw(value_error,"please use a step size of at least 1");
						}
						if (m_mode == inactive) {
							piranha_throw(value_error,"cannot calculate the limit of a power series expansion "
								"if no degree limit has been set");
						}
						if (s.empty()) {
							piranha_throw(value_error,"cannot calculate the limit of the power series expansion of "
								"an empty power series");
						}
						// min_degree will be either total or partial, depending on the mode.
						mp_rational min_degree(0);
						switch (m_mode) {
							case deg:
								min_degree = s.min_degree();
								break;
							case p_deg:
								min_degree = s.base_partial_min_degree(psyms2pos(m_psyms,args_tuple));
								break;
							case inactive:
								piranha_assert(false);
						}
						if (min_degree <= 0) {
							piranha_throw(value_error,"cannot calculate the limit of a power series expansion if the (partial) minimum degree "
								"of the series is negative or zero");
						}
						if (m_degree_limit < 0) {
							piranha_throw(value_error,"cannot calculate the limit of a power series expansion "
								"if the minimum degree limit is negative");
						}
						// (mp_rational(limit) / min_degree - start) / step_size + 1;
						mp_rational tmp(m_degree_limit);
						tmp /= min_degree;
						tmp -= start;
						tmp /= step_size;
						tmp += 1;
						if (tmp > 0) {
							// Take the floor of tmp and convert to integer.
							const int retval = (tmp.get_num() / tmp.get_den()).to_int();
							if (tmp == retval) {
								// If tmp was an integer in the beginning, we want to take the number
								// immediately preceding it (or zero).
								return (retval > 0) ? (retval - 1) : 0;
							} else {
								// If tmp was not an integer, let's take the floor.
								return retval;
							}
						} else {
							__PDEBUG(std::cout << "Negative power series limit calculated, inserting 0 instead." << '\n');
							return 0;
						}
					}
					bool is_effective() const
					{
						switch (m_mode) {
							case inactive:
								return false;
							case deg:
								return true;
							case p_deg:
								// In case of partial degree truncator is effective only if the position tuple
								// contains some elements.
								return (m_pos_tuple.template get<expo_args_pos>().size() > 0);
						}
						piranha_assert(false);
						return false;
					}
				protected:
					void init()
					{
						switch (m_mode) {
							case deg:
								std::sort(m_multiplier.m_terms1.begin(), m_multiplier.m_terms1.end(), min_degree_comparison());
								std::sort(m_multiplier.m_terms2.begin(), m_multiplier.m_terms2.end(), min_degree_comparison());
								break;
							case p_deg:
								// We need to do the sorting only if the position tuple
								// contains some elements.
								if (m_pos_tuple.template get<expo_args_pos>().size() > 0) {
									std::sort(m_multiplier.m_terms1.begin(), m_multiplier.m_terms1.end(),
										partial_min_degree_comparison(m_pos_tuple));
									std::sort(m_multiplier.m_terms2.begin(), m_multiplier.m_terms2.end(),
										partial_min_degree_comparison(m_pos_tuple));
								}
								break;
							case inactive:
								;
						}
					}
				private:
					Multiplier		&m_multiplier;
					pos_tuple_type		m_pos_tuple;
			};
			static void set(const int &);
			static void set(const mp_rational &);
			static void set(const vector_psym &, const int &);
			static void set(const vector_psym &, const mp_rational &);
			static void unset();
			static void print(std::ostream &stream = std::cout);
		private:
			static mp_rational	m_degree_limit;
			static vector_psym	m_psyms;
			static mode		m_mode;
	};

	typedef toolbox<degree_truncator_> degree_truncator;
}

#endif
