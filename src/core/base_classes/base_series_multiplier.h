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
			friend class base_insert_multiplication_result;
		protected:
			// Alias for term type of first input series and return value series.
			typedef typename Series1::term_type term_type1;
			// Alias for term type of second input series.
			typedef typename Series2::term_type term_type2;
			p_static_check((boost::is_same<typename term_type1::key_type, typename term_type2::key_type>::value),
				"Key type mismatch in base multiplier.");
			class key_revlex_comparison
			{
				public:
					template <class Term>
					bool operator()(const Term *t1, const Term *t2) const
					{
						return (t1->m_key.revlex_comparison(t2->m_key));
					}
			};
			template <class Term>
			struct ptr_cf_extractor
			{
				const typename Term::cf_type &operator()(const Term *t) {
					return t->m_cf;
				}
			};
			/// Compute block size for multiplication.
			/**
			 * Resulting block size depends on the cache memory available, and will always be in the [16,8192] range. N is
			 * the size of the base storage unit type used to store the result of the multiplication (e.g., coefficient type in vector coded, coded hash table
			 * term in the sparse hashed multiplication, etc.).
			 */
			template <std::size_t N>
			static std::size_t compute_block_size()
			{
				// NOTE: this function is used typically considering only the output storage requirements, since storage of input series
				//       will be a small fraction of storage for the output series.
				p_static_check(N > 0,"");
				const std::size_t shift = boost::numeric_cast<std::size_t>(
					std::log(std::max<double>(16.,std::sqrt((settings::cache_size * 1024) / N))) / std::log(2.) - 1
				);
				return (std::size_t(2) << std::min<std::size_t>(std::size_t(12),shift));
			}
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
		public:
			base_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
				m_s1(s1), m_s2(s2), m_args_tuple(args_tuple), m_size1(m_s1.length()),
				m_size2(m_s2.length()), m_retval(retval)
			{
				piranha_assert(m_size1 > 0 && m_size2 > 0);
			}
			// Plain multiplication.
			void perform_plain_multiplication()
			{
				perform_plain_threaded_multiplication();
			}
		private:
			template <class GenericTruncator>
			struct plain_functor {
				typedef typename term_type1::multiplication_result mult_res;
				plain_functor(mult_res &res,const term_type1 **t1, const term_type2 **t2, const GenericTruncator &trunc,
					Series1 &retval, const ArgsTuple &args_tuple):m_res(res),m_t1(t1),m_t2(t2),
					m_trunc(trunc),m_retval(retval),m_args_tuple(args_tuple)
				{}
				bool operator()(const std::size_t &i, const std::size_t &j)
				{
					if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
						return false;
					}
					term_type1::multiply(*m_t1[i], *m_t2[j], m_res, m_args_tuple);
					insert_multiplication_result<mult_res>::run(m_res, m_retval, m_args_tuple);
					return true;
				}
				mult_res		&m_res;
				const term_type1	**m_t1;
				const term_type2	**m_t2;
				const GenericTruncator	&m_trunc;
				Series1			&m_retval;
				const ArgsTuple		&m_args_tuple;
			};
			struct plain_worker {
				plain_worker(base_series_multiplier &mult, Series1 &retval):
					m_mult(mult),m_retval(retval),m_terms1(mult.m_terms1)
				{}
				plain_worker(base_series_multiplier &mult, Series1 &retval,
					std::vector<std::vector<term_type1 const *> > &split1, const std::size_t &idx):
					m_mult(mult),m_retval(retval),m_terms1(split1[idx])
				{}
				void operator()()
				{
					// Build the truncator.
					const typename Truncator::template get_type<Series1,Series2,ArgsTuple> trunc(m_terms1,m_mult.m_terms2,m_mult.m_args_tuple);
					// Use the selected truncator only if it really truncates, otherwise use the
					// null truncator.
					if (trunc.is_effective()) {
						plain_implementation(trunc);
					} else {
						plain_implementation(
							null_truncator::template get_type<Series1,Series2,ArgsTuple>(
							m_terms1,m_mult.m_terms2,m_mult.m_args_tuple
							)
						);
					}
				}
				template <class GenericTruncator>
				void plain_implementation(const GenericTruncator &trunc)
				{
					typedef typename term_type1::multiplication_result mult_res;
					mult_res res;
					const std::size_t size1 = m_terms1.size(), size2 = m_mult.m_size2;
					const term_type1 **t1 = &m_terms1[0];
					const term_type2 **t2 = &m_mult.m_terms2[0];
					plain_functor<GenericTruncator> pf(res,t1,t2,trunc,m_retval,m_mult.m_args_tuple);
					const std::size_t block_size = compute_block_size
						<boost::tuples::length<mult_res>::value * sizeof(term_type1)>();
					blocked_multiplication(block_size,size1,size2,pf);
				}
				base_series_multiplier		&m_mult;
				Series1				&m_retval;
				std::vector<term_type1 const *>	&m_terms1;
			};
			// Threaded multiplication.
			void perform_plain_threaded_multiplication()
			{
				// Effective number of threads to use. If the two series are small, we want to use one single thread.
				// NOTE: here the number 2500 is a kind of rule-of thumb. Basically multiplications of series < 50 elements
				// will use just one thread.
				if (double(m_size1) * double(m_size2) < 2500) {
					plain_worker w(*derived_cast,m_retval);
					w();
				} else {
					// If size1 is less than the number of desired threads,
					// use size1 as number of threads.
					const std::size_t n = std::min(settings::get_nthread(),m_size1);
					std::vector<std::vector<term_type1 const *> > split1(n);
					// m is the number of terms per thread for regular blocks.
					const std::size_t m = m_size1 / n;
					// Iterate up to n - 1 because that's the number up to which we can divide series1 into
					// regular blocks.
					for (std::size_t i = 0; i < n - 1; ++i) {
						split1[i].insert(split1[i].end(),m_terms1.begin() + i * m, m_terms1.begin() + (i + 1) * m);
					}
					// Last iteration.
					split1[n - 1].insert(split1[n - 1].end(),m_terms1.begin() + (n - 1) * m, m_terms1.end());
					boost::thread_group tg;
					std::vector<Series1> retvals(n,Series1());
					for (std::size_t i = 0; i < n; ++i) {
						tg.create_thread(plain_worker(*derived_cast,retvals[i],split1,i));
					}
					tg.join_all();
					// Take the retvals and insert them into final retval.
					for (std::size_t i = 0; i < n; ++i) {
						m_retval.insert_range(retvals[i].begin(),retvals[i].end(),m_args_tuple);
					}
				}
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
			// Vectors of pointers to the input terms.
			std::vector<term_type1 const *>			m_terms1;
			std::vector<term_type2 const *>			m_terms2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
