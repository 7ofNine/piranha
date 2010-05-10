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
#include <boost/lambda/lambda.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath> // For std::ceil.
#include <cstddef>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../psym.h"
#include "../settings.h" // For debug messages.

namespace piranha { namespace truncators {
	/// Truncator based on the minium degree of the series.
	class __PIRANHA_VISIBLE degree
	{
			enum mode {
				deg,
				p_deg,
				inactive
			};
			template <int ExpoTermPos>
			struct order_comparison
			{
				template <class Term>
				bool operator()(const Term *t1, const Term *t2) const
				{
					typedef typename Term::template component<ExpoTermPos>::type::degree_type degree_type;
					const degree_type md1(t1->template get<ExpoTermPos>().order()), md2(t2->template get<ExpoTermPos>().order());
					return md1 < md2;
				}
			};
			template <int ExpoTermPos, class PosTuple>
			struct partial_order_comparison
			{
				partial_order_comparison(const PosTuple &pos_tuple):m_pos_tuple(pos_tuple) {}
				template <class Term>
				bool operator()(const Term *t1, const Term *t2) const
				{
					typedef typename Term::template component<ExpoTermPos>::type::degree_type degree_type;
					const degree_type md1(t1->template get<ExpoTermPos>().partial_order(m_pos_tuple)),
						md2(t2->template get<ExpoTermPos>().partial_order(m_pos_tuple));
					return md1 < md2;
				}
				const PosTuple &m_pos_tuple;
			};
		public:
			template <class Series1, class Series2, class ArgsTuple>
			class get_type
			{
					typedef typename ntuple<std::vector<std::pair<bool,std::size_t> >,
						boost::tuples::length<ArgsTuple>::value>::type
						pos_tuple_type;
					static const int expo_term_pos = Series1::expo_term_position;
					static const int expo_args_pos = Series1::expo_args_position;
				public:
					typedef typename Series1::term_type term_type1;
					typedef typename Series2::term_type term_type2;
					typedef get_type type;
					get_type(std::vector<term_type1 const *> &t1, std::vector<term_type2 const *> &t2, const ArgsTuple &args_tuple, bool initialise = true):
						m_t1(t1),m_t2(t2),m_args_tuple(args_tuple)
					{
						// Some static checks.
						p_static_check(Series1::expo_args_position == Series2::expo_args_position, "");
						p_static_check(Series1::expo_term_position == Series2::expo_term_position, "");
						// Convert psyms vector into position tuple only if we are truncating to partial degree.
						if (m_mode == p_deg) {
							m_pos_tuple = psyms2pos(m_psyms,m_args_tuple);
						}
						if (initialise) {
							init();
						}
					}
					bool skip(const term_type1 **t1, const term_type2 **t2) const
					{
						switch (m_mode) {
							case deg:
								return ((*t1)->template get<expo_term_pos>().order() +
									(*t2)->template get<expo_term_pos>().order() >=
									m_degree_limit);
							case p_deg:
								return ((*t1)->template get<expo_term_pos>().partial_order(m_pos_tuple) +
									(*t2)->template get<expo_term_pos>().partial_order(m_pos_tuple) >=
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
					template <class PowerSeries, class ArgsTuple2>
					static std::size_t power_series_iterations(const PowerSeries &s, const int &start, const int &step_size,
						const ArgsTuple2 &args_tuple)
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
						// order will be either total or partial, depending on the mode.
						mp_rational order(0);
						switch (m_mode) {
							case deg:
								order = s.order();
								break;
							case p_deg:
								order = s.base_partial_order(psyms2pos(m_psyms,args_tuple));
								break;
							case inactive:
								piranha_assert(false);
						}
						if (order <= 0) {
							piranha_throw(value_error,"cannot calculate the limit of a power series expansion if the (partial) minimum degree "
								"of the series is negative or zero");
						}
						if (m_degree_limit < 0) {
							piranha_throw(value_error,"cannot calculate the limit of a power series expansion "
								"if the minimum degree limit is negative");
						}
						// (mp_rational(limit) / order - start) / step_size + 1;
						mp_rational tmp(m_degree_limit);
						tmp /= order;
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
					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::term_type const *> get_sorted_pointer_vector(const Series &s, const ArgsTuple2 &args_tuple)
					{
						std::vector<typename Series::term_type const *> retval;
						std::transform(s.begin(),s.end(),std::insert_iterator<std::vector<typename Series::term_type const *> >(retval,retval.begin()),
							&boost::lambda::_1);
						switch (m_mode) {
							case deg:
								std::sort(retval.begin(),retval.end(),order_comparison<Series::expo_term_position>());
								break;
							case p_deg:
								{
								typedef typename ntuple<std::vector<std::pair<bool,std::size_t> >,
									boost::tuples::length<ArgsTuple2>::value>::type pos_tuple_type;
								const pos_tuple_type pos_tuple(psyms2pos(m_psyms,args_tuple));
								if (pos_tuple.template get<Series::expo_args_position>().size() > 0) {
									std::sort(retval.begin(), retval.end(),partial_order_comparison<Series::expo_term_position,pos_tuple_type>(pos_tuple));
								} else {
									piranha_throw(value_error,"cannot establish series ordering, partial degree truncator is not "
										"effective on this series");
								}
								}
								break;
							case inactive:
								piranha_throw(value_error,"cannot establish series ordering, degree truncator is not active");
						}
						return retval;
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
								std::sort(m_t1.begin(), m_t1.end(), order_comparison<expo_term_pos>());
								std::sort(m_t2.begin(), m_t2.end(), order_comparison<expo_term_pos>());
								break;
							case p_deg:
								// We need to do the sorting only if the position tuple
								// contains some elements.
								if (m_pos_tuple.template get<expo_args_pos>().size() > 0) {
									std::sort(m_t1.begin(), m_t1.end(),
										partial_order_comparison<expo_term_pos,pos_tuple_type>(m_pos_tuple));
									std::sort(m_t2.begin(), m_t2.end(),
										partial_order_comparison<expo_term_pos,pos_tuple_type>(m_pos_tuple));
								}
								break;
							case inactive:
								;
						}
					}
				private:
					std::vector<term_type1 const *>	&m_t1;
					std::vector<term_type2 const *>	&m_t2;
					const ArgsTuple			&m_args_tuple;
					pos_tuple_type			m_pos_tuple;
			};
			static void set(const int &);
			static void set(const mp_rational &);
			static void set(const std::string &, const int &);
			static void set(const std::string &, const mp_rational &);
			static void set(const std::vector<std::string> &, const int &);
			static void set(const std::vector<std::string> &, const mp_rational &);
			static void unset();
			static void print(std::ostream &stream = std::cout);
		private:
			static mp_rational	m_degree_limit;
			static vector_psym	m_psyms;
			static mode		m_mode;
	};
} }

#endif
