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

#ifndef PIRANHA_EXPO_TRUNCATOR_H
#define PIRANHA_EXPO_TRUNCATOR_H

#include <algorithm> // To calculate min degree and for sorting.
#include <cmath> // For std::ceil.
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"
#include "../psym.h"
#include "../settings.h" // For debug messages.

namespace piranha
{
	/// Base exponent truncator.
	/**
	 * Internally it stores a static list of piranha::psym_p + int pairs representing the limits
	 * for the exponents of symbols. This class can be used by those classes that need the
	 * exponent limits (such as truncators for polynomial multiplication, classes that perform power series
	 * expansions, etc.).
	 */
	class __PIRANHA_VISIBLE base_expo_truncator
	{
		protected:
			typedef std::vector<std::pair<psym_p, max_fast_int> > container_type;
			typedef container_type::iterator iterator;
		public:
			static void print(std::ostream &stream = std::cout);
			static void set(const std::string &name, const max_fast_int &n) {
				generic_limit(name, n);
			}
			static void set(const psym &p, const max_fast_int &n) {
				generic_limit(p, n);
			}
			static void unset();
			static void unset(const std::string &name) {
				generic_unset(name);
			}
			static void unset(const psym &p) {
				generic_unset(p);
			}
			// Limit of a power series development of a power series.
			template <class PowerSeries, class ArgsTuple>
			static size_t power_series_limit(const PowerSeries &s, const ArgsTuple &args_tuple,
											 const int &start = 0, const int &step_size = 1) {
				p_assert(start >= 0 and step_size >= 1);
				if (s.empty()) {
					throw unsuitable("Cannot calculate the limit of the power series expansion of an empty power series.");
				}
				// Let's calculate the minimum exponents of s for which the truncator defines a limit.
				const std::pair<std::vector<max_fast_int>, std::vector<max_fast_int> >
				mle(min_limited_exponents(s, args_tuple));
				const size_t size = mle.first.size();
				if (size == 0) {
					throw not_existing("Cannot calculate the limit of a power series expansion when there are no exponent limits "
									   "set for the arguments of the series.");
				}
				for (size_t i = 0; i < size; ++i) {
					// If the minimum degree of the symbol whose exponent we want to limit is zero or less, we are
					// screwed, since the degree of the symbol would not increase at every step of the expansion
					// and we would end up in an infinite loop.
					if (mle.first[i] <= 0) {
						throw unsuitable("Cannot calculate the limit of a power series expansion if one of the limited exponents "
										 "of the series has negative or zero minimum value.");
					}
					if (mle.second[i] < 0) {
						throw unsuitable("Cannot calculate the limit of a power series expansion if there are negative "
										 "exponent limits for this series.");
					}
				}
				std::vector<size_t> power_limits(size);
				for (size_t i = 0; i < size; ++i) {
					const double tmp((static_cast<double>(mle.second[i]) / mle.first[i] - start) / static_cast<double>(step_size));
					if (tmp >= 0) {
						power_limits[i] = (size_t)std::ceil(tmp);
					} else {
						__PDEBUG(std::cout << "Negative power series limit calculated, inserting 0 instead." << '\n');
						power_limits[i] = 0;
					}
				}
				// The biggest limit is the one that defines the length of the power series expansion.
				size_t retval = *(std::max_element(power_limits.begin(), power_limits.end()));
				__PDEBUG(std::cout << "Calculated limit for power series of power series: " << retval << '\n');
				return retval;
			}
		protected:
			/// Transform the list of psymbol limits into a list of positions - limits given a piranha::vector_psym_p.
			static std::vector<std::pair<size_t, max_fast_int> > positions_limits(const vector_psym_p &v) {
				std::vector<std::pair<size_t, max_fast_int> > retval;
				const size_t limits_size = m_expo_limits.size(),
										   args_size = v.size();
				for (size_t i = 0; i < limits_size; ++i) {
					for (size_t j = 0; j < args_size; ++j) {
						if (m_expo_limits[i].first == v[j]) {
							retval.push_back(std::pair<size_t, max_fast_int>(j, m_expo_limits[i].second));
							// We can break out, there should not be duplicates inside the arguments list.
							break;
						}
					}
				}
				return retval;
			}
		private:
			template <class Argument>
			static void generic_limit(const Argument &arg, const max_fast_int &n) {
				psym_p tmp(psyms::get_pointer(arg));
				iterator it = find_argument(tmp);
				if (it == m_expo_limits.end()) {
					m_expo_limits.push_back(std::pair<psym_p, max_fast_int>(tmp, n));
				} else {
					it->second = n;
				}
			}
			template <class Argument>
			static void generic_unset(const Argument &p) {
				psym_p tmp(psyms::get_pointer(p));
				iterator it = find_argument(tmp);
				if (it == m_expo_limits.end()) {
					throw(not_existing(std::string("Symbol ") + "\"" + tmp->name() + "\" does not have an exponent limit set."));
				} else {
					m_expo_limits.erase(it);
				}
			}
			// Returns minimum_exponent - limit pairs.
			template <class PowerSeries, class ArgsTuple>
			static std::pair<std::vector<max_fast_int>, std::vector<max_fast_int> >
			min_limited_exponents(const PowerSeries &s, const ArgsTuple &args_tuple) {
				// First let's find out the minimum exponents of the input series.
				const std::vector<max_fast_int> min_expos(s.min_exponents(args_tuple));
				// Secondly, find out the position of the limited exponents.
				const std::vector<std::pair<size_t, max_fast_int> >
				pos_lim(positions_limits(args_tuple.template get<PowerSeries::expo_args_position>()));
				const size_t size = pos_lim.size();
				std::pair<std::vector<max_fast_int>, std::vector<max_fast_int> > retval;
				retval.first.resize(size);
				retval.second.resize(size);
				for (size_t i = 0; i < size; ++i) {
					p_assert(pos_lim[i].first < min_expos.size());
					retval.first[i] = min_expos[pos_lim[i].first];
					retval.second[i] = pos_lim[i].second;
				}
				return retval;
			}
			static iterator find_argument(const psym_p &);
		protected:
			static container_type m_expo_limits;
	};

	/// Truncators for polynomials based on the exponent of one or more variables.
	class expo_truncator
	{
		public:
			template <class Multiplier>
			class get_type: public base_expo_truncator
			{
					static const int expo_term_pos = Multiplier::series_type1::expo_term_position;
					static const int expo_args_pos = Multiplier::series_type1::expo_args_position;
					class single_expo_comparison
					{
							typedef typename Multiplier::args_tuple_type args_tuple_type;
						public:
							single_expo_comparison(const size_t &pos, const args_tuple_type &args_tuple):
									m_pos(pos), m_args_tuple(args_tuple) {}
							template <class Term>
							bool operator()(const Term &t1, const Term &t2) const {
								return t1.template get<expo_term_pos>().min_expo_of(m_pos, m_args_tuple) <
									   t2.template get<expo_term_pos>().min_expo_of(m_pos, m_args_tuple);
							}
						private:
							const size_t			&m_pos;
							const args_tuple_type	&m_args_tuple;
					};
				public:
					typedef get_type type;
					get_type(Multiplier &m):
							m_multiplier(m),
							m_positions(base_expo_truncator::positions_limits(
											m_multiplier.m_args_tuple.template
											get<expo_args_pos>())),
							m_size(m_positions.size()) {
						// Some static checks.
						p_static_check(Multiplier::series_type1::expo_args_position ==
									   Multiplier::series_type2::expo_args_position, "");
						p_static_check(Multiplier::series_type1::expo_term_position ==
									   Multiplier::series_type2::expo_term_position, "");
						// If just one exponent is limited, sort series according the the min degree of such exponent.
						if (m_positions.size() == 1) {
							std::sort(m_multiplier.m_terms1.begin(), m_multiplier.m_terms1.end(),
									  single_expo_comparison(m_positions[0].first, m_multiplier.m_args_tuple));
							std::sort(m_multiplier.m_terms2.begin(), m_multiplier.m_terms2.end(),
									  single_expo_comparison(m_positions[0].first, m_multiplier.m_args_tuple));
						}
					}
					bool accept(const max_fast_int &n) const {
						if (m_size <= 1) {
							return true;
						} else {
							m_multiplier.decode(m_multiplier.m_tmp_key, n);
							return m_multiplier.m_tmp_key.test_expo_limits(m_positions, m_multiplier.m_args_tuple);
						}
					}
					template <class Term>
					bool accept(const Term &t) const {
						if (m_size <= 1) {
							return true;
						} else {
							return t.template get<expo_term_pos>().test_expo_limits(
									   m_positions, m_multiplier.m_args_tuple
								   );
						}
					}
					template <class Term1, class Term2>
					bool skip(const Term1 &t1, const Term2 &t2) const {
						if (m_size == 1 &&
								(t1.template get<expo_term_pos>().min_expo_of(m_positions[0].first, m_multiplier.m_args_tuple) +
								 t2.template get<expo_term_pos>().min_expo_of(m_positions[0].first, m_multiplier.m_args_tuple) >
								 m_positions[0].second)) {
							return true;
						}
						return false;
					}
					bool is_effective() const {
						return !m_positions.empty();
					}
				private:
					Multiplier											&m_multiplier;
					const std::vector<std::pair<size_t, max_fast_int> >	m_positions;
					const size_t										m_size;
			};
	};
}

#endif
