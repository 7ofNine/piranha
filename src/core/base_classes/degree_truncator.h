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
#include <cmath> // For std::ceil.
#include <iostream>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"
#include "../settings.h" // For debug messages.

namespace piranha
{
	/// Base degree truncator.
	class __PIRANHA_VISIBLE base_degree_truncator
	{
		public:
			static void set(const max_fast_int &n) {
				m_degree_limit = n;
				m_effective = true;
			}
			static void clear();
			static void print(std::ostream &stream = std::cout);
			// Limit of a power series development of a power series according to its minimum degree.
			template <class PowerSeries, class ArgsTuple>
			static size_t power_series_limit(const PowerSeries &s, const ArgsTuple &,
											 const int &start = 0, const int &step_size = 1) {
				p_assert(start >= 0 and step_size >= 1);
				if (!m_effective) {
					throw unsuitable("Cannot calculate the limit of a power series expansion if no degree limit has been set.");
				}
				if (s.empty()) {
					throw unsuitable("Cannot calculate the limit of the power series expansion of an empty power series.");
				}
				const max_fast_int min_degree(s.min_degree());
				if (min_degree <= 0) {
						throw unsuitable("Cannot calculate the limit of a power series expansion if the minimum degree "
										 "of the series is negative or zero.");
				}
				if (m_degree_limit < 0) {
					throw unsuitable("Cannot calculate the limit of a power series expansion if the minimum degree limit "
										"is negative.");
				}
				const double tmp((static_cast<double>(m_degree_limit) / min_degree - start) / static_cast<double>(step_size));
				if (tmp >= 0) {
					return static_cast<size_t>(std::ceil(tmp));
				} else {
					__PDEBUG(std::cout << "Negative power series limit calculated, inserting 0 instead." << '\n');
					return 0;
				}
			}
		protected:
			static max_fast_int	m_degree_limit;
			static bool			m_effective;
	};

	/// Truncator based on the minium degree of the series.
	class degree_truncator
	{
		public:
			template <class Multiplier>
			class get_type: public base_degree_truncator
			{
					static const int expo_term_pos = Multiplier::series_type1::expo_term_position;
					static const int expo_args_pos = Multiplier::series_type1::expo_args_position;
					class min_degree_comparison
					{
						public:
							template <class Term>
							bool operator()(const Term &t1, const Term &t2) const {
								return t1.template get<expo_term_pos>().min_degree() <
									   t2.template get<expo_term_pos>().min_degree();
							}
					};
				public:
					typedef get_type type;
					get_type(Multiplier &m): m_multiplier(m) {
						// Some static checks.
						p_static_check(Multiplier::series_type1::expo_args_position ==
									   Multiplier::series_type2::expo_args_position, "");
						p_static_check(Multiplier::series_type1::expo_term_position ==
									   Multiplier::series_type2::expo_term_position, "");
						// Sort series according to the minimum degree.
						std::sort(m_multiplier.m_terms1.begin(), m_multiplier.m_terms1.end(),min_degree_comparison());
						std::sort(m_multiplier.m_terms2.begin(), m_multiplier.m_terms2.end(),min_degree_comparison());
					}
					template <class T>
					bool accept(const T &) const {
						return true;
					}
					template <class Term1, class Term2>
					bool skip(const Term1 &t1, const Term2 &t2) const {
						return (t1.template get<expo_term_pos>().min_degree() + t2.template get<expo_term_pos>().min_degree() >
								 m_degree_limit);
					}
				private:
					Multiplier	&m_multiplier;
			};
	};
}

#endif
