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

#ifndef PIRANHA_CODED_SERIES_MULTIPLIER_H
#define PIRANHA_CODED_SERIES_MULTIPLIER_H

#include <boost/integer_traits.hpp> // For integer limits.
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <vector>

#include "../p_assert.h"
#include "../settings.h" // For debug messages.

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Toolbox for coded series multiplication.
	/**
	 * Intended to be inherited together with piranha::base_series_multiplier. It adds common methods for coded
	 * series multiplication.
	 */
	template <class Derived>
	class coded_series_multiplier
	{
			typedef boost::integer_traits<max_fast_int> traits;
		public:
			// Decode code n into key.
			template <class Key>
			void decode(Key &key, const max_fast_int &n) const {
				key.decode(n, m_coding_vector, m_h_min, m_fast_res_min_max, derived_const_cast->m_args_tuple);
			}
		protected:
			coded_series_multiplier():
					m_cr_is_viable(false),
					m_size(derived_const_cast->m_args_tuple.template get<Derived::key_type::position>().size()),
					m_min_max1(m_size), m_min_max2(m_size), m_res_min_max(m_size), m_fast_res_min_max(m_size),
					// Coding vector is larger to accomodate extra element at the end (used during decodification).
					m_coding_vector(m_size + 1) {
				m_ckeys1.reserve(derived_const_cast->m_size1);
				m_ckeys2.reserve(derived_const_cast->m_size2);
			}
			void find_input_min_max() {
				size_t i1 = 0, i2 = 0;
				// Fill first minmax vector. This works because at this point we are sure both series have
				// at least one term. Assert it, just to make sure.
				p_assert(derived_const_cast->m_size1 > 0 && derived_const_cast->m_size2 > 0);
				derived_const_cast->m_terms1[i1].m_key.upload_ints_to(m_min_max1);
				derived_const_cast->m_terms2[i2].m_key.upload_ints_to(m_min_max2);
				// Move to the second terms and cycle on all remaining terms.
				++i1;
				++i2;
				for (; i1 < derived_const_cast->m_size1; ++i1) {
					derived_cast->m_terms1[i1].m_key.test_min_max_ints(m_min_max1);
				}
				for (; i2 < derived_const_cast->m_size2; ++i2) {
					derived_cast->m_terms2[i2].m_key.test_min_max_ints(m_min_max2);
				}
				__PDEBUG(std::cout << "Limits are:\n";
				for (size_t i = 0; i < m_min_max1.size(); ++i) {
				std::cout << m_min_max1[i].first << ',' << m_min_max1[i].second << '\n';
				}
				std::cout << "and:\n";
				for (size_t i = 0; i < m_min_max2.size(); ++i) {
				std::cout << m_min_max2[i].first << ',' << m_min_max2[i].second << '\n';
				})
			}
			void determine_viability() {
				// We must do the computations with arbitrary integers to avoid exceeding range.
				mpz_class hmin(0), hmax(0), ck(1);
				for (size_t i = 0; i < m_size; ++i) {
					hmin += ck * m_res_min_max[i].first;
					hmax += ck * m_res_min_max[i].second;
					// Assign also the coding vector, so we avoid doing it later.
					m_coding_vector[i] = ck.get_si();
					ck *= (m_res_min_max[i].second - m_res_min_max[i].first + 1);
				}
				__PDEBUG(std::cout << "hmax-hmin=" << hmax - hmin << '\n');
				// We want to fill an extra slot of the coding vector (wrt to the nominal size,
				// corresponding to the arguments number for the key). This is handy for decodification.
				m_coding_vector[m_size] = ck.get_si();
				p_assert(ck > 0);
				// Determine viability by checking that ck and the minimum/maximum values for the codes
				// respect the fast integer boundaries.
				if (ck < traits::const_max && hmin > traits::const_min && hmin < traits::const_max &&
						hmax > traits::const_min && hmax < traits::const_max) {
					m_cr_is_viable = true;
					m_h_min = hmin.get_si();
					m_h_max = hmax.get_si();
					// Downcast minimum and maximum result values to fast integers.
					for (size_t i = 0; i < m_size; ++i) {
						if (m_res_min_max[i].first < traits::const_min || m_res_min_max[i].first > traits::const_max ||
								m_res_min_max[i].second < traits::const_min || m_res_min_max[i].second > traits::const_max) {
							std::cout << "Warning: results of series multiplication cross " <<
									  "fast integer limits. Expect errors." << std::endl;
						}
						m_fast_res_min_max[i].first = m_res_min_max[i].first.get_si();
						m_fast_res_min_max[i].second = m_res_min_max[i].second.get_si();
					}
					__PDEBUG(std::cout << "Coding vector: ";
					for (size_t i = 0; i < m_size; ++i) {
					std::cout << m_coding_vector[i] << '\t';
					}
					std::cout << "+\t" << m_coding_vector[m_size] << '\n';)
				}
			}
			/// Code keys.
			void code_keys() {
				for (size_t i = 0; i < derived_const_cast->m_size1; ++i) {
					m_ckeys1.push_back(derived_const_cast->m_terms1[i].m_key.code(m_coding_vector,
									   derived_const_cast->m_args_tuple));
				}
				for (size_t i = 0; i < derived_const_cast->m_size2; ++i) {
					m_ckeys2.push_back(derived_const_cast->m_terms2[i].m_key.code(m_coding_vector,
									   derived_const_cast->m_args_tuple));
				}
			}
		protected:
			// Is coded representation viable?
			bool													m_cr_is_viable;
			// Size of the coding vector, min_max vectors, etc.
			const size_t											m_size;
			// Vectors of minimum and maximum value pairs for the series being multiplied.
			std::vector<std::pair<max_fast_int, max_fast_int> >		m_min_max1;
			std::vector<std::pair<max_fast_int, max_fast_int> >		m_min_max2;
			// Vector of minimum and maximum value pairs for the resulting series.
			// GMP is used to avoid trespassing the range limits of max_fast_int.
			std::vector<std::pair<mpz_class, mpz_class> >			m_res_min_max;
			// Version of the above downcast to fast integer type.
			std::vector<std::pair<max_fast_int, max_fast_int> >		m_fast_res_min_max;
			// Coding vector.
			std::vector<max_fast_int>								m_coding_vector;
			// Mininum and maximum values of codes.
			max_fast_int											m_h_min;
			max_fast_int											m_h_max;
			// Coded keys.
			std::vector<max_fast_int>								m_ckeys1;
			std::vector<max_fast_int>								m_ckeys2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
