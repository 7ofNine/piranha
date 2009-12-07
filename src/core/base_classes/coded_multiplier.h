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

#ifndef PIRANHA_CODED_MULTIPLIER_H
#define PIRANHA_CODED_MULTIPLIER_H

#include <boost/integer_traits.hpp>
#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <exception>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"
#include "../ntuple.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Series, int N>
	struct minmax_tuple {
		p_static_check(N > 0,"");
		typedef typename Series::term_type::key_type::value_type value_type;
		typedef boost::tuples::cons<std::vector<std::pair<value_type,value_type> >,
			typename minmax_tuple<typename Series::term_type::cf_type,N - 1>::type> type;
	};

	template <class Series>
	struct minmax_tuple<Series,0> {
		typedef boost::tuples::null_type type;
	};

	template <class Derived, class Series1, class Series2, class ArgsTuple, class OpsTuple>
	class coded_series_multiplier {
			p_static_check(boost::tuples::length<ArgsTuple>::value == boost::tuples::length<OpsTuple>::value);
			typedef typename minmax_tuple<Series::> minmax_type
			typedef typename ntuple<std::vector<mp_integer>,boost::tuples::length<ArgsTuple>::value>::type tuple_type;
		protected:
			
			tuple_type 
	};


	/// Toolbox for coded series multiplication.
	/**
	 * Intended to be inherited together with piranha::base_series_multiplier. It adds common methods for
	 * series multiplication through Kronecker codification. Requirements:
	 * - input series must have at least one term.
	 */
	template <class Derived, class Key, bool SubtractKeys>
	class coded_series_multiplier
	{
		protected:
			/// Key value type.
			typedef typename Key::value_type value_type;
			/// Key size type.
			typedef typename Key::size_type size_type;
			/// Default constructor.
			/**
			 * Initialises viable flag to false and sets up data members for future use.
			 */
			// NOTE: check the order of inheritance here?
			coded_multiplier():m_cr_is_viable(false),m_size(derived_const_cast->m_args_tuple.template get<Key::position>().size()),
				m_min_max1(m_size),m_min_max2(m_size),m_res_min_max(m_size),m_coding_vector(),m_h_min(0),m_h_max(0),m_h_tot(0),
				m_ckeys1(derived_const_cast->m_size1),m_ckeys2(derived_const_cast->m_size2),m_density1(0.),m_density2(0.)
			{
				if (m_size == boost::integer_traits<size_type>::const_max) {
					piranha_throw(std::overflow_error,"overflow in key's width during codification");
				}
				m_coding_vector.resize(m_size + size_type(1));
			}
			void find_input_min_max()
			{
				std::size_t i1 = 0, i2 = 0;
				const std::size_t size1 = derived_const_cast->m_size1, size2 = derived_const_cast->m_size2;
				// Fill first minmax vector. This works because at this point we are sure both series have
				// at least one term. Assert it, just to make sure.
				piranha_assert(size1 > 0 && size2 > 0);
				piranha_assert(derived_const_cast->m_terms1[i1]->m_key.size() == m_size);
				for (size_type i = 0; i < derived_const_cast->m_terms1[i1]->m_key.size(); ++i) {
					m_min_max1[i].first = m_min_max1[i].second =
						derived_const_cast->m_terms1[i1]->m_key[i];
				}
				piranha_assert(derived_const_cast->m_terms2[i2]->m_key.size() <= m_size);
				for (size_type i = 0; i < derived_const_cast->m_terms2[i2]->m_key.size(); ++i) {
					m_min_max2[i].first = m_min_max2[i].second =
						derived_const_cast->m_terms2[i2]->m_key[i];
				}
				// Move to the second terms and cycle on all remaining terms.
				++i1;
				++i2;
				for (; i1 < size1; ++i1) {
					derived_cast->m_terms1[i1]->m_key.test_min_max_ints(m_min_max1);
				}
				for (; i2 < size2; ++i2) {
					derived_cast->m_terms2[i2]->m_key.test_min_max_ints(m_min_max2);
				}
				__PDEBUG(std::cout << "Limits are:\n";
				for (std::size_t i = 0; i < m_min_max1.size(); ++i) {
				std::cout << m_min_max1[i].first << ',' << m_min_max1[i].second << '\n';
				}
				std::cout << "and:\n";
				for (std::size_t i = 0; i < m_min_max2.size(); ++i) {
				std::cout << m_min_max2[i].first << ',' << m_min_max2[i].second << '\n';
				})
			}
		protected:
			/// Is coded representation viable?
			bool							m_cr_is_viable;
			/// Size of the coding vector, min_max vectors, etc.
			const size_type						m_size;
			/// Vectors of minimum and maximum value pairs for the first series.
			std::vector<std::pair<value_type, value_type> >		m_min_max1;
			/// Vectors of minimum and maximum value pairs for the second series.
			std::vector<std::pair<value_type, value_type> >		m_min_max2;
			/// Vector of minimum and maximum value pairs for the resulting series.
			std::vector<std::pair<value_type, value_type> >		m_res_min_max;
			/// Coding vector.
			std::vector<max_fast_int>				m_coding_vector;
			/// Mininum code.
			max_fast_int						m_h_min;
			/// Maximum code.
			max_fast_int						m_h_max;
			/// Total number of codes.
			max_fast_int						m_h_tot;
			/// Coded keys for the first series.
			std::vector<max_fast_int>				m_ckeys1;
			/// Coded keys for the second series.
			std::vector<max_fast_int>				m_ckeys2;
			/// Density of the first series.
			double							m_density1;
			/// Density of the second series.
			double							m_density2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
