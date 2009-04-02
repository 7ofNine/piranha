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

#ifndef PIRANHA_NORM_TRUNCATOR_H
#define PIRANHA_NORM_TRUNCATOR_H

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../config.h"
#include "../exceptions.h"
#include "../none.h"
#include "../p_assert.h"
#include "../settings.h"

namespace piranha
{
	class __PIRANHA_VISIBLE base_norm_truncator
	{
		public:
			static void set(const int &n) {
				if (n <= 0) {
					throw(unsuitable("Please insert a positive integer."));
				} else {
					m_truncation_level = std::pow(10., -n);
				}
				m_truncation_power = n;
			}
			static void print(std::ostream &stream = std::cout);
			// Returns the length of a development in powers of x that satisfies the condition that the
			// magnitude of the last term of the expansion with respect to x's magnitudes is m_truncation_level
			// times smaller.
			template <class T, class ArgsTuple>
			static size_t power_series_iterations(const T &x, const int &start, const int &step_size,
				const ArgsTuple &args_tuple) {
				// NOTE: share this check in some kind of base truncator class?
				if (step_size < 1) {
					throw unsuitable("Please use a step size of at least 1.");
				}
				if (m_truncation_power == 0) {
					throw unsuitable("No value set for norm-based truncation, "
						"cannot calculate limit of power series expansion.");
				}
				const double norm = x.norm_(args_tuple);
				p_assert(norm >= 0);
				if (norm >= 1) {
					throw unsuitable("The norm of the argument of the power series expansion is >= 1: "
						"the norm truncator is unable to give an estimate of the power series limit.");
				}
				// Let's prevent log10(0) below.
				if (norm == 0) {
					throw unsuitable("Unable to find a limit for the power series expansion of a series whose norm "
						"is zero.");
				}
				int retval = static_cast<int>(std::ceil(static_cast<double>(
						static_cast<int>(std::ceil(std::log10(m_truncation_level) /
						std::log10(norm) + 1 - start)) / step_size))) + 1;
				// This could be negative if starting power is big enough. In this case return 0.
				if (retval >= 0) {
					return retval;
				} else {
					return 0;
				}
			}
			static void unset();
			static bool is_effective() {
				return m_truncation_power != 0;
			}
		protected:
			static int	m_truncation_power;
			static double	m_truncation_level;
	};

	/// Norm-based truncator.
	class norm_truncator
	{
		public:
			template <class Multiplier>
			class get_type: public base_norm_truncator
			{
					template <class ArgsTuple>
					class norm_comparison
					{
						public:
							norm_comparison(const ArgsTuple &args_tuple): m_args_tuple(args_tuple) {}
							template <class Term>
							bool operator()(const Term *t1, const Term *t2) const {
								return (t1->m_cf.norm_(m_args_tuple) * t1->m_key.norm_(m_args_tuple) >
										t2->m_cf.norm_(m_args_tuple) * t2->m_key.norm_(m_args_tuple));
							}
						private:
							const ArgsTuple	&m_args_tuple;
					};
				public:
					typedef get_type type;
					get_type(Multiplier &m, bool initialise = true):
							m_multiplier(m),
							m_delta_threshold(
								m.m_s1.norm_(m.m_args_tuple)*m.m_s2.norm_(m.m_args_tuple)*m_truncation_level /
								(2*m.m_s1.length()*m.m_s2.length())) {
						if (initialise) {
							init();
						}
					}
					template <class Result>
					bool accept(const Result &) const {
						return true;
					}
					template <class Term1, class Term2>
					bool skip(const Term1 &t1, const Term2 &t2) const {
						return (
								t1.m_cf.norm_(m_multiplier.m_args_tuple) *
								t1.m_key.norm_(m_multiplier.m_args_tuple) *
								t2.m_cf.norm_(m_multiplier.m_args_tuple) *
								t2.m_key.norm_(m_multiplier.m_args_tuple) / 2. <
								m_delta_threshold
						);
					}
				protected:
					void init() {
						if (is_effective()) {
							const norm_comparison<typename Multiplier::args_tuple_type> cmp(m_multiplier.m_args_tuple);
							std::sort(m_multiplier.m_terms1.begin(), m_multiplier.m_terms1.end(), cmp);
							std::sort(m_multiplier.m_terms2.begin(), m_multiplier.m_terms2.end(), cmp);
						}
					}
				private:
					Multiplier    &m_multiplier;
					const double  m_delta_threshold;
			};
	};
}

#endif
