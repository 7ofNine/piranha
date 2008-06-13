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

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/type_traits/is_same.hpp> // For key type detection.
#include <vector>

#include "../p_assert.h"
#include "../proxies.h"
#include "../settings.h"
#include "base_series_multiplier_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Base series multiplier.
	/**
	 * This class is meant to be extended to build specific multipliers.
	 */
	template <class Series1, class Series2, class ArgsTuple, template <class> class Truncator, class Derived>
	class base_series_multiplier
	{
			friend struct base_insert_multiplication_result;
			friend class Truncator<Derived>;
		protected:
			// Alias for term type of first input series and return value series.
			typedef typename Series1::term_type term_type1;
			// Alias for term type of second input series.
			typedef typename Series2::term_type term_type2;
			// Alias for the truncator type.
			typedef Truncator<Derived> truncator_type;
		private:
			typedef boost::multi_index_container
			<
			term_type1,
			boost::multi_index::indexed_by
			<
			boost::multi_index::hashed_unique<boost::multi_index::identity<term_type1> >
			>
			>
			mult_set;
			BOOST_STATIC_ASSERT((boost::is_same<typename term_type1::key_type, typename term_type2::key_type>::value));
			typedef cf_mult_proxy<typename term_type1::cf_type> cf_proxy_type1;
			typedef cf_mult_proxy<typename term_type2::cf_type> cf_proxy_type2;
			typedef key_mult_proxy<typename term_type1::key_type> key_proxy_type;
		public:
			base_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
					m_s1(s1), m_s2(s2), m_args_tuple(args_tuple), m_size1(m_s1.template nth_index<0>().size()),
					m_size2(m_s2.template nth_index<0>().size()), m_retval(retval), m_terms1(m_size1), m_terms2(m_size2),
					m_trunc(*derived_cast) {
				// Set proper load factor for hash set.
				m_set.max_load_factor(settings::load_factor());
				// Cache pointers the terms of the series into vectors.
				cache_series_terms(m_s1, m_terms1);
				cache_series_terms(m_s2, m_terms2);
			}
		protected:
			/// Plain multiplication.
			void perform_plain_multiplication() {
				plain_multiplication();
				plain_insert_result_into_retval();
			}
		private:
			// Perform plain multiplication.
			void plain_multiplication() {
				typedef typename term_type1::multiplication_result mult_res;
				mult_res res;
				for (size_t i = 0; i < m_size1; ++i) {
					for (size_t j = 0; j < m_size2; ++j) {
						if (m_trunc.skip(m_terms1[i], m_terms2[j])) {
							break;
						}
						term_type1::multiply(m_terms1[i], m_terms2[j], res, m_args_tuple);
						insert_multiplication_result<mult_res>::run(res, *this);
					}
				}
			}
			template <class Series>
			void cache_series_terms(const Series &s,
									std::vector<typename Series::term_type const *> &terms) {
				typedef typename Series::const_sorted_iterator const_sorted_iterator;
				const const_sorted_iterator it_f = s.template nth_index<0>().end();
				size_t i = 0;
				for (const_sorted_iterator it = s.template nth_index<0>().begin(); it != it_f; ++it) {
					terms[i] = &(*it);
					++i;
				}
			}
			// After the multiplication has been performed and the result stored in the temporary hash table,
			// fetch the terms from there and put them into retval.
			void plain_insert_result_into_retval() {
				typedef typename mult_set::const_iterator hash_iterator;
				typedef typename Series1::sorted_iterator sorted_iterator;
				term_type1 term;
				sorted_iterator it_hint = m_retval.template nth_index<0>().end();
				const hash_iterator it_f = m_set.end();
				for (hash_iterator it = m_set.begin(); it != it_f; ++it) {
					term.m_cf = it->m_cf;
					term.m_key = it->m_key;
					it_hint = m_retval.template insert<false, true>(term, it_hint, m_args_tuple);
				}
			}
		protected:
			// References to the series.
			const Series1                 	&m_s1;
			const Series2                 	&m_s2;
			// Reference to the arguments tuple.
			const ArgsTuple               	&m_args_tuple;
			// Sizes of the series.
			const size_t                  	m_size1;
			const size_t                  	m_size2;
			// Reference to the result.
			Series1                       	&m_retval;
			// Vectors of pointers to the input terms.
			std::vector<term_type1 const *>	m_terms1;
			std::vector<term_type2 const *>	m_terms2;
			// Container to store the result of the multiplications.
			mult_set                      	m_set;
			// Truncator. This must be the last one defined because it will take *this
			// as parameter for construction.
			truncator_type                	m_trunc;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
