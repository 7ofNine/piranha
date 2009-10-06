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
#include <boost/thread/thread.hpp>
#include <boost/type_traits/is_same.hpp> // For key type detection.
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <cstddef>
#include <vector>

#include "../base_classes/null_truncator.h"
#include "../config.h"
#include "../runtime.h"
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
			class key_revlex_comparison
			{
				public:
					template <class Term>
					bool operator()(const Term *t1, const Term *t2) const
					{
						return (t1->m_key.revlex_comparison(t2->m_key));
					}
			};
			template <class Cf>
			struct cf_from_term {
				template <class Term>
				static const Cf &get(const Term *t)
				{
					return t->m_cf;
				}
			};
			template <class Cf>
			struct cf_direct {
				static const Cf &get(const Cf &c)
				{
					return c;
				}
			};
			template <class Functor>
			static void blocked_multiplication(const std::size_t &block_size, const std::size_t &size1, const std::size_t &size2, Functor &m)
			{
				piranha_assert(block_size > 0);
				const std::size_t nblocks1 = size1 / block_size, nblocks2 = size2 / block_size;
				for (std::size_t n1 = 0; n1 < nblocks1; ++n1) {
					const std::size_t i_start = n1 * block_size, i_end = i_start + block_size;
					// regulars1 * regulars2
					for (std::size_t n2 = 0; n2 < nblocks2; ++n2) {
						const std::size_t j_start = n2 * block_size, j_end = j_start + block_size;
						for (std::size_t i = i_start; i < i_end; ++i) {
							for (std::size_t j = j_start; j < j_end; ++j) {
								if (!m(i,j)) {
									break;
								}
							}
						}
					}
					// regulars1 * rem2
					for (std::size_t i = i_start; i < i_end; ++i) {
						for (std::size_t j = nblocks2 * block_size; j < size2; ++j) {
							if (!m(i,j)) {
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
							if (!m(i,j)) {
								break;
							}
						}
					}
				}
				// rem1 * rem2.
				for (std::size_t i = nblocks1 * block_size; i < size1; ++i) {
					for (std::size_t j = nblocks2 * block_size; j < size2; ++j) {
						if (!m(i,j)) {
							break;
						}
					}
				}
			}
			template <class GenericTruncator>
			struct plain_functor {
				typedef typename term_type1::multiplication_result mult_res;
				plain_functor(mult_res &res,const term_type1 **t1, const term_type2 **t2, const GenericTruncator &trunc,
					Series1 &retval, const ArgsTuple &args_tuple):m_res(res),m_t1(t1),m_t2(t2),
					m_trunc(trunc),m_retval(retval),m_args_tuple(args_tuple) {}
				bool operator()(const std::size_t &i, const std::size_t &j)
				{
					if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
						return false;
					}
					term_type1::multiply(*m_t1[i], *m_t2[j], m_res, m_args_tuple);
					insert_multiplication_result<mult_res>::run(m_res, m_retval, m_trunc, m_args_tuple);
					return true;
				}
				mult_res		&m_res;
				const term_type1	**m_t1;
				const term_type2	**m_t2;
				const GenericTruncator	&m_trunc;
				Series1			&m_retval;
				const ArgsTuple		&m_args_tuple;
			};
			template <class Multiplier>
			struct plain_worker {
				plain_worker(Multiplier &mult, Series1 &retval, const std::size_t &idx):
					m_mult(mult),m_retval(retval),m_idx(idx) {}
				void operator()()
				{
					// Build the truncator.
					const typename Truncator::template get_type<Series1,Series2,ArgsTuple> trunc(m_mult.m_split1[m_idx],m_mult.m_split2[m_idx],m_mult.m_args_tuple);
					// Use the selected truncator only if it really truncates, otherwise use the
					// null truncator.
					if (trunc.is_effective()) {
						plain_implementation(trunc);
					} else {
						plain_implementation(
							null_truncator::template get_type<Series1,Series2,ArgsTuple>(
							m_mult.m_split1[m_idx],m_mult.m_split2[m_idx],m_mult.m_args_tuple
							)
						);
					}
				}
				template <class GenericTruncator>
				void plain_implementation(const GenericTruncator &trunc)
				{
					typedef typename term_type1::multiplication_result mult_res;
					mult_res res;
					const std::size_t size1 = m_mult.m_split1[m_idx].size(), size2 = m_mult.m_size2;
					const term_type1 **t1 = &m_mult.m_split1[m_idx][0];
					const term_type2 **t2 = &m_mult.m_split2[m_idx][0];
					plain_functor<GenericTruncator> pf(res,t1,t2,trunc,m_retval,m_mult.m_args_tuple);
					const std::size_t block_size = 2 << (
							(std::size_t)log2(std::max(16.,std::sqrt((settings::cache_size * 1024) /
							((sizeof(term_type1) + sizeof(term_type2) + boost::tuples::length<mult_res>::value * sizeof(term_type1))
							* runtime::get_n_cur_threads())))) - 1);
std::cout << "Block size: " << block_size << '\n';
					blocked_multiplication(block_size,size1,size2,pf);
				}
				Multiplier			&m_mult;
				Series1				&m_retval;
				const std::size_t		m_idx;
			};
		private:
			p_static_check((boost::is_same<typename term_type1::key_type, typename term_type2::key_type>::value),
				"Key type mismatch in base multiplier.");
		public:
			base_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
				m_s1(s1), m_s2(s2), m_args_tuple(args_tuple), m_size1(m_s1.length()),
				m_size2(m_s2.length()), m_retval(retval),
				m_terms1(utils::cache_terms_pointers(s1)),m_terms2(utils::cache_terms_pointers(s2)),
				m_split1()
			{
				piranha_assert(m_size1 > 0);
				// Effective number of threads to use. If the two series are small, we want to use one single thread.
				// TODO: testing to see which number should go here. Maybe test with small poly multiplication and see how many of them we
				// can do per second with the best possible scenario (double coefficients, integer exponents, vector coded) and compare to how
				// many threads we can generate per second with little overhead.
				// TODO: think about dropping m_terms* altogether and using only m_split*.
				std::size_t n;
				if (double(m_size1) * double(m_size2) < 400) {
					n = 1;
				} else {
					// If size1 is less than the number of desired threads,
					// use size1 as number of threads.
					n = std::min(settings::get_nthread(),m_size1);
				}
				piranha_assert(n > 0);
				m_split1.reserve(n);
				// m is the number of terms per thread for regular blocks.
				const std::size_t m = m_size1 / n;
				// Iterate up to n - 1 because that's the number up to which we can divide series1 into
				// regular blocks.
				for (std::size_t i = 0;i < n - 1; ++i) {
					m_split1.push_back(std::vector<term_type1 const *>(m_terms1.begin() + i * m,m_terms1.begin() + (i + 1) * m));
					m_split2.push_back(m_terms2);
					std::sort(m_split1[i].begin(),m_split1[i].end(),key_revlex_comparison());
					std::sort(m_split2[i].begin(),m_split2[i].end(),key_revlex_comparison());
				}
				// Last iteration.
				m_split1.push_back(std::vector<term_type1 const *>(m_terms1.begin() + (n - 1) * m,m_terms1.end()));
				m_split2.push_back(m_terms2);
			}
			// Threaded multiplication.
			template <class Worker>
			void perform_threaded_multiplication()
			{
				const std::size_t n = m_split1.size();
				piranha_assert(n > 0);
				if (n == 1) {
					Worker w(*derived_cast,m_retval,0);
					w();
				} else {
std::cout << "Going threaded\n";
					boost::thread_group tg;
					std::vector<Series1> retvals(n,Series1());
					{
					runtime::register_threads r(n - 1);
					for (std::size_t i = 0; i < n; ++i) {
						tg.create_thread(Worker(*derived_cast,retvals[i],i));
					}
std::cout << "joining\n";
					tg.join_all();
					}
std::cout << "joined\n";
					// Take the retvals and insert them into final retval.
					for (std::size_t i = 0; i < n; ++i) {
						m_retval.insert_range(retvals[i].begin(),retvals[i].end(),m_args_tuple);
					}
std::cout << "inserted\n";
				}
			}
			// Plain multiplication.
			void perform_plain_multiplication()
			{
				perform_threaded_multiplication<plain_worker<base_series_multiplier> >();
			}
		public:
			// References to the series.
			const Series1					&m_s1;
			const Series2					&m_s2;
			// Reference to the arguments tuple.
			const ArgsTuple					&m_args_tuple;
			// Sizes of the series.
			const std::size_t				m_size1;
			const std::size_t				m_size2;
			// Reference to the result.
			Series1						&m_retval;
			// Vectors of pointers the input terms.
			std::vector<term_type1 const *>			m_terms1;
			std::vector<term_type2 const *>			m_terms2;
			// Vector resulting from splitting m_terms1 into chunks and copies of m_terms2 to be used in threads.
			std::vector<std::vector<term_type1 const *> >	m_split1;
			std::vector<std::vector<term_type2 const *> >	m_split2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
