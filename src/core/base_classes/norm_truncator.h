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

#include <cmath>
#include <iostream>
#include <string>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"

namespace piranha
{
	class __PIRANHA_VISIBLE base_norm_truncator
	{
		public:
			static void set(const int &n) {
				if (n < 0) {
					throw(unsuitable("Please insert a non-negative integer."));
				} else if (n == 0) {
					m_truncation_level = 0;
				} else {
					m_truncation_level = std::pow(10., -n);
				}
			}
			static void print(std::ostream &stream = std::cout);
			static std::string print_to_string();
			// Returns the length of a development in powers of x that satisfies the condition that the
			// magnitude of the last term of the expansion with respect to x's magnitudes is m_truncation_level
			// times smaller.
			template <class T, class ArgsTuple>
			static size_t power_series_limit(const T &x, const ArgsTuple &args_tuple,
											 const int &start = 0, const int &step_size = 1) {
				p_assert(step_size >= 1 and start >= 0);
				const double norm = x.norm(args_tuple);
				p_assert(norm >= 0);
				if (norm >= 1) {
					throw unsuitable("The norm of the argument of the power series expansion is >= 1: the expansion will diverge.");
				}
				max_fast_int retval = (max_fast_int)std::ceil((max_fast_int)std::ceil(std::log10(m_truncation_level)
									  / std::log10(norm) + 1 - start) / step_size);
				if (retval >= 0) {
					return retval;
				} else {
					return 0;
				}
			}
		private:
			static double m_truncation_level;
	};

	/// Norm-based truncator.
	template <class Multiplier>
	class norm_truncator: public base_norm_truncator
	{
		public:
			norm_truncator(Multiplier &m):
					m_multiplier(m),
					m_delta_threshold(
						m.m_s1.norm(m.m_args_tuple)*m.m_s2.norm(m.m_args_tuple)*m_truncation_level /
						(2*m.m_s1.template nth_index<0>().size()*m.m_s2.template nth_index<0>().size())) {}
			template <class Result>
			bool accept(const Result &) const {
				return true;
			}
			template <class Cf1, class Cf2, class Key>
			bool skip(const Cf1 &c1, const Key &, const Cf2 &c2, const Key &) const {
				return (c1.norm(m_multiplier.m_args_tuple) * c2.norm(m_multiplier.m_args_tuple) / 2 < m_delta_threshold);
			}
		private:
			Multiplier    &m_multiplier;
			const double  m_delta_threshold;
	};
}

#endif
