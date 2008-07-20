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

#include "../config.h"
#include "../p_assert.h"
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
	template <class Series1, class Series2, class ArgsTuple, class Truncator, class Derived>
	class base_series_multiplier
	{
			friend struct base_insert_multiplication_result;
			friend class Truncator::template get_type<Derived>::type;
		protected:
			// Alias for term type of first input series and return value series.
			typedef typename Series1::term_type term_type1;
			// Alias for term type of second input series.
			typedef typename Series2::term_type term_type2;
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
			p_static_check((boost::is_same<typename term_type1::key_type, typename term_type2::key_type>::value),
				"Key type mismatch in base multiplier.");
			typedef typename term_type1::template rebind < typename term_type1::cf_type::proxy::type,
			typename term_type1::key_type::proxy::type >::type term_proxy_type1;
			typedef typename term_type2::template rebind < typename term_type2::cf_type::proxy::type,
			typename term_type2::key_type::proxy::type >::type term_proxy_type2;
		public:
			base_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
					m_s1(s1), m_s2(s2), m_args_tuple(args_tuple), m_size1(m_s1.length()),
					m_size2(m_s2.length()), m_retval(retval),
					m_terms1(m_s1.cache_proxies()),m_terms2(m_s2.cache_proxies()) {
				// Set proper load factor for hash set.
				m_set.max_load_factor(settings::load_factor());
			}
		protected:
			/// Plain multiplication.
			template <class GenericTruncator>
			void perform_plain_multiplication(const GenericTruncator &trunc) {
				plain_multiplication(trunc);
				plain_insert_result_into_retval();
			}
		private:
			// Perform plain multiplication.
			template <class GenericTruncator>
			void plain_multiplication(const GenericTruncator &trunc) {
				typedef typename term_type1::multiplication_result mult_res;
				mult_res res;
				for (size_t i = 0; i < m_size1; ++i) {
					for (size_t j = 0; j < m_size2; ++j) {
						if (trunc.skip(m_terms1[i], m_terms2[j])) {
							break;
						}
						term_type1::multiply(m_terms1[i], m_terms2[j], res, m_args_tuple);
						insert_multiplication_result<mult_res>::run(res, *this, trunc);
					}
				}
			}
			// After the multiplication has been performed and the result stored in the temporary hash table,
			// fetch the terms from there and put them into retval.
			void plain_insert_result_into_retval() {
				typedef typename mult_set::const_iterator hash_iterator;
				term_type1 term;
				const hash_iterator it_f = m_set.end();
				for (hash_iterator it = m_set.begin(); it != it_f; ++it) {
					term.m_cf = it->m_cf;
					term.m_key = it->m_key;
					m_retval.template insert<false, true>(term, m_args_tuple);
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
			// Vectors of proxies for the input terms.
			std::vector<term_proxy_type1>	m_terms1;
			std::vector<term_proxy_type2>	m_terms2;
			// Container to store the result of the multiplications.
			mult_set                      	m_set;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
