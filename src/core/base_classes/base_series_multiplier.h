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

#ifndef PIRANHA_BASE_SERIES_MULTIPLIER_H
#define PIRANHA_BASE_SERIES_MULTIPLIER_H

#include <algorithm>
#include <boost/type_traits/is_same.hpp> // For key type detection.
#include <cstddef>
#include <utility>
#include <vector>

#include "../config.h"
#include "../settings.h"
#include "../utils.h"
#include "base_series_multiplier_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Base series multiplier.
	/**
	 * This class is meant to be extended to build specific multipliers.
	 */
	template <class Series1, class Series2, class ArgsTuple, class Truncator, class Derived>
	class base_series_multiplier
	{
			friend class base_insert_multiplication_result;
		protected:
			// Alias for term type of first input series and return value series.
			typedef typename Series1::term_type term_type1;
			// Alias for term type of second input series.
			typedef typename Series2::term_type term_type2;
			// Alias for vector of ranges used during parallel multiplication.
			typedef std::vector<std::pair<typename std::vector<term_type1 const *>::const_iterator,typename std::vector<term_type1 const *>::const_iterator> >
				vector_ranges;
			class key_revlex_comparison
			{
				public:
					template <class Term>
					bool operator()(const Term *t1, const Term *t2) const {
						return (t1->m_key.revlex_comparison(t2->m_key));
					}
			};
			template <class Cf>
			struct cf_from_term {
				template <class Term>
				static const Cf &get(const Term *t) {
					return t->m_cf;
				}
			};
			template <class Cf>
			struct cf_direct {
				static const Cf &get(const Cf &c) {
					return c;
				}
			};
			template <std::size_t block_size, template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
				class Term1, class Term2, class Ckey, class Trunc, class Result, class Multiplier>
			static void blocked_multiplication(const std::size_t &size1, const std::size_t &size2,
				const TermOrCf1 *tc1, const TermOrCf2 *tc2, const Term1 **t1, const Term2 **t2,
				const Ckey *ck1, const Ckey *ck2, const Trunc &trunc, Result *res, Multiplier &m,
				const ArgsTuple &args_tuple)
			{
				p_static_check(block_size > 0, "Invalid block size for cache-blocking.");
				const std::size_t nblocks1 = size1 / block_size, nblocks2 = size2 / block_size;
				for (std::size_t n1 = 0; n1 < nblocks1; ++n1) {
					const std::size_t i_start = n1 * block_size, i_end = i_start + block_size;
					// regulars1 * regulars2
					for (std::size_t n2 = 0; n2 < nblocks2; ++n2) {
						const std::size_t j_start = n2 * block_size, j_end = j_start + block_size;
						for (std::size_t i = i_start; i < i_end; ++i) {
							for (std::size_t j = j_start; j < j_end; ++j) {
								if (!m.template run<CfGetter>(i,j,tc1,tc2,t1,t2,ck1,ck2,trunc,res,args_tuple)) {
									break;
								}
							}
						}
					}
					// regulars1 * rem2
					for (std::size_t i = i_start; i < i_end; ++i) {
						for (std::size_t j = nblocks2 * block_size; j < size2; ++j) {
							if (!m.template run<CfGetter>(i,j,tc1,tc2,t1,t2,ck1,ck2,trunc,res,args_tuple)) {
								break;
							}
						}
					}
				}
				// rem1 * regulars2
				for (std::size_t n2 = 0; n2 < nblocks2; ++n2) {
					const std::size_t j_start = n2 * block_size, j_end = j_start + block_size;
					for (std::size_t i = nblocks1 * block_size; i < size1; ++i) {
						for (std::size_t j = j_start; j < j_end; ++j) {
							if (!m.template run<CfGetter>(i,j,tc1,tc2,t1,t2,ck1,ck2,trunc,res,args_tuple)) {
								break;
							}
						}
					}
				}
				// rem1 * rem2.
				for (std::size_t i = nblocks1 * block_size; i < size1; ++i) {
					for (std::size_t j = nblocks2 * block_size; j < size2; ++j) {
						if (!m.template run<CfGetter>(i,j,tc1,tc2,t1,t2,ck1,ck2,trunc,res,args_tuple)) {
							break;
						}
					}
				}
			}
		private:
			p_static_check((boost::is_same<typename term_type1::key_type, typename term_type2::key_type>::value),
				"Key type mismatch in base multiplier.");
		public:
			base_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
				m_s1(s1), m_s2(s2), m_args_tuple(args_tuple), m_size1(m_s1.length()),
				m_size2(m_s2.length()), m_retval(retval),
				m_terms1(utils::cache_terms_pointers(s1)),m_terms2(utils::cache_terms_pointers(s2))
			{
				piranha_assert(m_size1 > 0);
				// Effective number of threads to use. If size1 is less than the number of desired threads,
				// use size1 as number of threads.
				const std::size_t n = std::min(settings::get_nthread(),m_size1);
				piranha_assert(n > 0);
				m_ranges.reserve(n);
				// m is the number of terms per thread for homogeneous blocks.
				const std::size_t m = m_size1 / n;
				// Iterate up to n - 1 because that's the number up to which we can divide series1 into
				// homogeneous blocks.
				for (std::size_t i = 0;i < n - 1; ++i) {
					m_ranges.push_back(std::make_pair(m_terms1.begin() + i * m,m_terms1.begin() + (i + 1) * m));
				}
				// Last iteration.
				m_ranges.push_back(std::make_pair(m_terms1.begin() + (n - 1) * m,m_terms1.end()));
			}
			// Perform plain multiplication.
			template <class GenericTruncator>
			void perform_plain_multiplication(const GenericTruncator &trunc)
			{
				typedef typename term_type1::multiplication_result mult_res;
				mult_res res;
				for (std::size_t i = 0; i < m_size1; ++i) {
					for (std::size_t j = 0; j < m_size2; ++j) {
						if (trunc.skip(&m_terms1[i], &m_terms2[j])) {
							break;
						}
						term_type1::multiply(*m_terms1[i], *m_terms2[j], res, m_args_tuple);
						insert_multiplication_result<mult_res>::run(res, m_retval, trunc, m_args_tuple);
					}
				}
			}
		public:
			// References to the series.
			const Series1			&m_s1;
			const Series2			&m_s2;
			// Reference to the arguments tuple.
			const ArgsTuple			&m_args_tuple;
			// Sizes of the series.
			const std::size_t		m_size1;
			const std::size_t		m_size2;
			// Reference to the result.
			Series1				&m_retval;
			// Vectors of pointers the input terms.
			std::vector<term_type1 const *>	m_terms1;
			std::vector<term_type2 const *>	m_terms2;
			// Vectors of ranges
			vector_ranges			m_ranges;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
