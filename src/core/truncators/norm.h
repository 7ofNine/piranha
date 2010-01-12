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
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "../base_classes/toolbox.h"
#include "../config.h"
#include "../exceptions.h"
#include "../settings.h"
#include "../utils.h"

namespace piranha
{
	struct norm_;

	/// Norm-based truncator.
	template <>
	class __PIRANHA_VISIBLE toolbox<norm_>
	{
			template <class ArgsTuple>
			class norm_comparison
			{
				public:
					norm_comparison(const ArgsTuple &args_tuple): m_args_tuple(args_tuple) {}
					template <class Term>
					bool operator()(const Term *t1, const Term *t2) const
					{
						const double n1(t1->m_cf.norm(m_args_tuple) * t1->m_key.norm(m_args_tuple)), n2(t2->m_cf.norm(m_args_tuple) * t2->m_key.norm(m_args_tuple));
						if (n1 != n2) {
							return (n1 > n2);
						} else {
							// NOTICE: the idea is that for leading terms with equal
							// norm we choose the ones that have
							// unity key vector, so that we increase the chance of
							// being able to perform the expansion.
							if (t1->m_key.is_unity()) {
								return true;
							} else if (t2->m_key.is_unity()) {
								return false;
							}
							return (t1->m_key < t2->m_key);
						}

					}
				private:
					const ArgsTuple	&m_args_tuple;
			};
		public:
			template <class Series1, class Series2, class ArgsTuple>
			class get_type
			{
				public:
					typedef typename Series1::term_type term_type1;
					typedef typename Series2::term_type term_type2;
					typedef get_type type;
					get_type(std::vector<term_type1 const *> &terms1, std::vector<term_type2 const *> &terms2, const ArgsTuple &args_tuple, bool initialise = true):
						m_terms1(terms1),m_terms2(terms2),m_args_tuple(args_tuple),m_t1(0),m_t2(0)
					{
						if (initialise) {
							init();
						}
					}
					bool skip(const term_type1 **t1, const term_type2 **t2) const
					{
						return (m_t1 && t1 >= m_t1) || (m_t2 && t2 >= m_t2);
					}
					bool is_effective() const
					{
						return m_truncation_level != 0;
					}
					// Returns the length of a development in powers of x that satisfies the condition that the
					// magnitude of the last term of the expansion with respect to x's magnitudes is m_truncation_level
					// times smaller.
					template <class T, class ArgsTuple2>
					static std::size_t power_series_iterations(const T &x, const int &start, const int &step_size,
						const ArgsTuple2 &args_tuple)
					{
						// NOTE: share this check in some kind of base truncator class?
						if (step_size < 1) {
							piranha_throw(value_error,
								"please use a step size of at least 1");
						}
						if (m_truncation_level == 0) {
							piranha_throw(value_error,"no value set for norm-based truncation, "
								"cannot calculate limit of power series expansion");
						}
						const double norm = x.base_norm(args_tuple);
						piranha_assert(norm >= 0);
						if (norm >= 1) {
							piranha_throw(value_error,"the norm of the argument of the power series expansion is >= 1: "
								"the norm truncator is unable to give an estimate of the power series limit");
						}
						// Let's prevent log10(0) below.
						if (norm == 0) {
							piranha_throw(value_error,"unable to find a limit for the power series expansion of a series whose norm "
								"is zero");
						}
						const int retval = boost::numeric_cast<int>(std::ceil(
							(std::log(m_truncation_level) / std::log(norm) + step_size - start) / step_size + 1
						));
						// This could be negative if starting power is big enough. In this case return 0.
						return (retval >= 0) ? retval : 0;
					}
					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::term_type const *> get_sorted_pointer_vector(const Series &s, const ArgsTuple2 &args_tuple)
					{
						std::vector<typename Series::term_type const *> retval;
						utils::cache_terms_pointers(s,retval);
						if (m_truncation_level == 0) {
							piranha_throw(value_error,"cannot establish series ordering, norm truncator is not active");
						}
						const norm_comparison<ArgsTuple2> cmp(args_tuple);
						std::sort(retval.begin(), retval.end(), cmp);
						return retval;
					}
				protected:
					void init()
					{
						if (is_effective()) {
							piranha_assert(m_terms1.size() >= 1 && m_terms2.size() >= 1);
							const norm_comparison<ArgsTuple> cmp(m_args_tuple);
							std::sort(m_terms1.begin(), m_terms1.end(), cmp);
							std::sort(m_terms2.begin(), m_terms2.end(), cmp);
							const double norm1 = calculate_norm(m_terms1), norm2 = calculate_norm(m_terms2);
							// If one of the norms is zero, then we don't want to do any multiplication.
							const term_type1 **final1 = &(*(m_terms1.begin()));
							const term_type2 **final2 = &(*(m_terms2.begin()));
							if (norm1 == 0 || norm2 == 0) {
								m_t1 = final1;
								m_t2 = final2;
								return;
							}
							const double
								limit1 = -norm2 + std::sqrt(norm2 * norm2 + m_truncation_level * norm2 / norm1),
								limit2 = -norm1 + std::sqrt(norm1 * norm1 + m_truncation_level * norm1 / norm2);
							double delta1 = 0, delta2 = 0;
							// NOTE: the issue here is that it may happen that limit is greater than norm, in which case all the calculations
							// done for maximimzing the number of terms truncated do not apply. In such case for the time being we just turn
							// off truncation, but it should be probably dealt with differently.
							while (limit1 < norm1) {
								const term_type1 **tmp = (m_t1) ? m_t1 - 1 : &(*(m_terms1.end() - 1));
								delta1 += (*tmp)->m_cf.norm(m_args_tuple) * (*tmp)->m_key.norm(m_args_tuple);
								// If, by going to the next term, we exceed the delta1, leave m_t1 where it is and break out.
								if (delta1 >= limit1) {
									break;
								}
								// Assign new limit.
								m_t1 = tmp;
								// If we reached the top, break out.
								if (tmp == final1) {
									break;
								}
							}
							while (limit2 < norm2) {
								const term_type2 **tmp = (m_t2) ? m_t2 - 1 : &(*(m_terms2.end() - 1));
								delta2 += (*tmp)->m_cf.norm(m_args_tuple) * (*tmp)->m_key.norm(m_args_tuple);
								// If, by going to the next term, we exceed the delta2, leave m_t2 where it is and break out.
								if (delta2 >= limit2) {
									break;
								}
								// Assign new limit.
								m_t2 = tmp;
								// If we reached the top, break out.
								if (tmp == final2) {
									break;
								}
							}
						}
					}
				private:
					template <class Term>
					double calculate_norm(const std::vector<Term const *> &v) const
					{
						double retval = 0;
						const std::size_t end = v.size();
						for (std::size_t i = 0; i < end; ++i) {
							retval += v[i]->m_cf.norm(m_args_tuple) * v[i]->m_key.norm(m_args_tuple);
						}
						return retval;
					}
					std::vector<term_type1 const *>	&m_terms1;
					std::vector<term_type2 const *>	&m_terms2;
					const ArgsTuple			&m_args_tuple;
					term_type1 const		**m_t1;
					term_type2 const		**m_t2;
			};
			// Shared portion.
			static void set(const double &x)
			{
				if (x <= 0) {
					piranha_throw(value_error,"please insert a positive number");
				}
				m_truncation_level = x;
			}
			static void print(std::ostream &stream = std::cout);
			static void unset();
		private:
			static double m_truncation_level;
	};

namespace truncators
{
	typedef toolbox<norm_> norm;
}
}

#endif
