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

#ifndef PIRANHA_POLYNOMIAL_MULTIPLIER_H
#define PIRANHA_POLYNOMIAL_MULTIPLIER_H

#include <algorithm> // For std::max.
#include <boost/bind.hpp>
#include <boost/integer_traits.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/thread/thread.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <cstddef>
#include <exception>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_multiplier.h"
#include "../coded_hash_table.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../settings.h" // For debug and cache size.
#include "../stats.h"

namespace piranha
{
	// Threaded vector multiplication.
	template <class Functor>
	static inline void threaded_vector_blocked_multiplication(const std::size_t &block_size, const std::size_t &size1, const std::size_t &size2, const std::size_t &thread_id,
		const std::size_t &thread_n, boost::barrier *b, Functor &m)
	{
		piranha_assert(block_size > 0 && thread_n > 0 && thread_id < thread_n);
		// Numerical limits check. We need an extra block size buffer at the end to make sure we are able to
		// represent all indices and sizes.
		piranha_assert(size1 < boost::integer_traits<std::size_t>::const_max - block_size && size2 < boost::integer_traits<std::size_t>::const_max - block_size);
		// Number of blocks of dimension block_size.
		const std::size_t nrblocks1 = size1 / block_size, nrblocks2 = size2 / block_size;
		// Size of the last "irregular" two blocks.
		const std::size_t ib_size1 = size1 % block_size, ib_size2 = size2 % block_size;
		// Total number of blocks.
		const std::size_t nblocks1 = nrblocks1 + (ib_size1 != 0), nblocks2 = nrblocks2 + (ib_size2 != 0);
		// Number of block iterations: for the first series we want to jump every n_thread blocks
		// (also maybe going past the series end), while for the second series we want to iterate
		// over all blocks plus allow also for thread_n - 1 "silent" iterations which are needed
		// to make sure that at the same time the threads are operating on ordered blocks (which
		// guarantees that we are writing in isolated areas of the output array).
		const std::size_t nbi1 = nblocks1 / thread_n + ((nblocks1 % thread_n) != 0),
			nbi2 = nblocks2 + (thread_n - 1);
		for (std::size_t n1 = 0; n1 < nbi1; ++n1) {
			const std::size_t i_start = (thread_id + n1 * thread_n) * block_size;
			// Here the block size will be zero if we are past the end of the first series, ib_size1 if this is the last block,
			// block_size if this is a standard block.
			const std::size_t cur_block_size1 = (i_start >= size1) ? 0 : ((i_start + block_size > size1) ? ib_size1 : block_size);
			const std::size_t i_end = i_start + cur_block_size1;
			for (std::size_t n2 = 0; n2 < nbi2; ++n2) {
				const std::size_t j_start = ((thread_id + n2) % nbi2) * block_size;
				const std::size_t cur_block_size2 = (j_start >= size2) ? 0 : ((j_start + block_size > size2) ? ib_size2 : block_size);
				const std::size_t j_end = j_start + cur_block_size2;
				// NOTE: here maybe we can put a preemptive check on j start/end so that if the inner
				// cycle is empty we skip this part altogether.
				for (std::size_t i = i_start; i < i_end; ++i) {
					for (std::size_t j = j_start; j < j_end; ++j) {
						m(i,j);
					}
				}
				// Synchronize this thread.
				b->wait();
			}
		}
	}

	/// Series multiplier specifically tuned for polynomials.
	/**
	 * This multiplier internally will use coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
	class polynomial_multiplier
	{
		public:
			template <class Series1, class Series2, class ArgsTuple, class Truncator>
			class get_type:
				public base_series_multiplier< Series1, Series2, ArgsTuple, Truncator,
					get_type<Series1, Series2, ArgsTuple, Truncator> >,
				public coded_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,boost::tuple<boost::true_type> >
			{
					typedef base_series_multiplier< Series1, Series2, ArgsTuple, Truncator,
						get_type<Series1, Series2, ArgsTuple, Truncator> > ancestor;
					typedef coded_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,boost::tuple<boost::true_type> > coded_ancestor;
					friend class coded_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,boost::tuple<boost::true_type> >;
					typedef typename ancestor::term_type1 term_type1;
					typedef typename ancestor::term_type2 term_type2;
					typedef typename term_type1::cf_type cf_type1;
					typedef typename term_type2::cf_type cf_type2;
				public:
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<Series1,Series2,ArgsTuple> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
						ancestor(s1, s2, retval, args_tuple) {}
					template <class GenericTruncator>
					struct vector_functor {
						vector_functor(std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
							std::vector<max_fast_int> &ck1, std::vector<max_fast_int> &ck2,
							std::vector<term_type1 const *> &t1, std::vector<term_type2 const *> &t2,
							const GenericTruncator &trunc, cf_type1 *vc_res, const ArgsTuple &args_tuple):
							m_tc1(tc1),m_tc2(tc2),m_ck1(ck1),m_ck2(ck2),m_t1(t1),m_t2(t2),m_trunc(trunc),
							m_vc_res(vc_res),m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							// Calculate index of the result.
							const max_fast_int res_index = m_ck1[i] + m_ck2[j];
							m_vc_res[res_index].addmul(m_tc1[i],m_tc2[j],m_args_tuple);
							return true;
						}
						static void blocks_setup(const std::size_t &cur_idx1_start, const std::size_t &block_size,
							std::vector<std::size_t> &idx_vector1, std::vector<std::size_t> &idx_vector2)
						{
							// TODO: needs sorting if cur_idx1_start is zero.
						}
						std::vector<cf_type1>		&m_tc1;
						std::vector<cf_type2>		&m_tc2;
						std::vector<max_fast_int>	&m_ck1;
						std::vector<max_fast_int>	&m_ck2;
						std::vector<term_type1 const *> &m_t1;
						std::vector<term_type2 const *> &m_t2;
						const GenericTruncator		&m_trunc;
						cf_type1			*m_vc_res;
						const ArgsTuple			&m_args_tuple;
					};
					template <class GenericTruncator>
					bool perform_vector_coded_multiplication(std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
						std::vector<term_type1 const *> &t1, std::vector<term_type2 const *> &t2, const GenericTruncator &trunc)
					{
						std::vector<cf_type1,std_counting_allocator<cf_type1> > vc;
						// Try to allocate the space for vector coded multiplication.
						// The +1 is needed because we need the number of possible codes between min and max.
						piranha_assert(boost::numeric::width(this->m_fast_h) + 1 >= 0);
						const std::size_t n_codes = boost::numeric_cast<std::size_t>(boost::numeric::width(this->m_fast_h) + 1);
						try {
							vc.resize(n_codes);
						} catch (const std::bad_alloc &) {
							__PDEBUG(std::cout << "Not enough physical memory available for vector coded.\n");
							return false;
						} catch (const memory_error &) {
							__PDEBUG(std::cout << "Memory limit reached for vector coded.\n");
							return false;
						}
						__PDEBUG(std::cout << "Going for vector coded polynomial multiplication\n");
						// Define the base pointers for storing the results of multiplication.
						// NOTE: even if here it seems like we are going to write outside allocated memory,
						//       the indices from the analysis of the coded series will prevent out-of-boundaries
						//       reads/writes. The thing works like this: we have ncodes slots allocated, so memory
						//       indices in [0,ncodes - 1]. But since we are doing arithmetics on shifted codes, the
						//       code range is [-h_min,ncodes - 1 - h_min]. So we shift the baseline memory location
						//       so that we can use shifted codes directly as indices.
						const std::size_t size1 = this->m_terms1.size(), size2 = this->m_terms2.size();
						piranha_assert(size1 && size2);
						const args_tuple_type &args_tuple = this->m_args_tuple;
						cf_type1 *vc_res =  &vc[0] - this->m_fast_h.lower();
						// Find out a suitable block size.
						const std::size_t block_size = this->template compute_block_size<sizeof(cf_type1)>();
						__PDEBUG(std::cout << "Block size: " << block_size << '\n');
// std::cout << "Block size: " << block_size << '\n';
						// Perform multiplication.
						typedef vector_functor<GenericTruncator> vf_type;
						vf_type vm(tc1,tc2,this->m_ckeys1,this->m_ckeys2a,t1,t2,trunc,vc_res,args_tuple);
						const std::size_t nthread = settings::get_nthread();
// const boost::posix_time::ptime time0 = boost::posix_time::microsec_clock::local_time();
						if (trunc.is_effective() || (this->m_terms1.size() * this->m_terms2.size()) <= 400 || nthread == 1) {
							stats::trace_stat("mult_st",std::size_t(0),boost::lambda::_1 + 1);
							this->blocked_multiplication(block_size,size1,size2,vm);
						} else {
// std::cout << "using " << nthread << " threads\n";
							stats::trace_stat("mult_mt",std::size_t(0),boost::lambda::_1 + 1);
							boost::thread_group tg;
							boost::barrier b(nthread);
							for (std::size_t i = 0; i < nthread; ++i) {
								tg.create_thread(boost::bind(threaded_vector_blocked_multiplication<vf_type>,block_size,size1,size2,i,nthread,&b,vm));
							}
							tg.join_all();
						}
// std::cout << "Elapsed time: " << (double)(boost::posix_time::microsec_clock::local_time() - time0).total_microseconds() / 1000 << '\n';
						__PDEBUG(std::cout << "Done multiplying\n");
						const max_fast_int i_f = this->m_fast_h.upper();
						// Decode and insert the results into return value.
						term_type1 tmp_term;
						for (max_fast_int i = this->m_fast_h.lower(); i <= i_f; ++i) {
							// Take a shortcut and check for ignorability of the coefficient here.
							// This way we avoid decodification, and all the series term insertion yadda-yadda.
							if (!vc_res[i].is_ignorable(args_tuple)) {
								this->decode(vc_res[i], i,tmp_term);
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done polynomial vector coded.\n");
						return true;
					}
					template <class Cf, class Ckey, class GenericTruncator, class HashSet>
					struct hash_functor {
						hash_functor(std::pair<Cf,Ckey> &cterm, std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
							std::vector<Ckey> &ck1, std::vector<Ckey> &ck2,
							std::vector<const term_type1 *> &t1, std::vector<const term_type2 *> &t2,
							const GenericTruncator &trunc, HashSet *cms, const ArgsTuple &args_tuple):
							m_cterm(cterm),m_tc1(tc1),m_tc2(tc2),m_ck1(ck1),m_ck2(ck2),m_t1(t1),m_t2(t2),
							m_trunc(trunc),m_cms(cms),m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							typedef typename HashSet::iterator c_iterator;
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							m_cterm.second = m_ck1[i];
							m_cterm.second += m_ck2[j];
							std::pair<bool,c_iterator> res = m_cms->find(m_cterm.second);
							if (res.first) {
								res.second->first.addmul(m_tc1[i],m_tc2[j],m_args_tuple);
							} else {
								// Assign to the temporary term the old cf (new_key is already assigned).
								m_cterm.first = m_tc1[i];
								// Multiply the old term by the second term.
								m_cterm.first.mult_by(m_tc2[j],m_args_tuple);
								m_cms->insert_new(m_cterm,res.second);
							}
							return true;
						}
						std::pair<Cf,Ckey>		&m_cterm;
						std::vector<cf_type1>		&m_tc1;
						std::vector<cf_type2>		&m_tc2;
						std::vector<Ckey>		&m_ck1;
						std::vector<Ckey>		&m_ck2;
						std::vector<const term_type1 *>	&m_t1;
						std::vector<const term_type2 *> &m_t2;
						const GenericTruncator		&m_trunc;
						HashSet				*m_cms;
						const ArgsTuple			&m_args_tuple;
					};
					template <class GenericTruncator>
					void perform_hash_coded_multiplication(std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
						std::vector<const term_type1 *> &t1, std::vector<const term_type2 *> &t2, const GenericTruncator &trunc)
					{
						typedef coded_hash_table<cf_type1, max_fast_int, std_counting_allocator<char> > csht;
						typedef typename csht::iterator c_iterator;
						stats::trace_stat("mult_st",std::size_t(0),boost::lambda::_1 + 1);
						// Let's find a sensible size hint.
						const std::size_t n_codes = boost::numeric_cast<std::size_t>(boost::numeric::width(this->m_fast_h) + 1);
						const std::size_t size_hint = static_cast<std::size_t>(
							std::max<double>(this->m_density1,this->m_density2) * n_codes);
						const std::size_t size1 = this->m_terms1.size(), size2 = this->m_terms2.size();
						piranha_assert(size1 && size2);
						const args_tuple_type &args_tuple = this->m_args_tuple;
						csht cms(size_hint);
						// Find out a suitable block size.
						const std::size_t block_size = this->template compute_block_size<sizeof(std::pair<cf_type1,max_fast_int>)>();
						__PDEBUG(std::cout << "Block size: " << block_size << '\n');
// std::cout << "Block size: " << block_size << '\n';
// const boost::posix_time::ptime time0 = boost::posix_time::microsec_clock::local_time();
						std::pair<cf_type1,max_fast_int> cterm;
						hash_functor<cf_type1,max_fast_int,GenericTruncator,csht>
							hm(cterm,tc1,tc2,this->m_ckeys1,this->m_ckeys2a,t1,t2,trunc,&cms,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,hm);
//std::cout << "Elapsed time: " << (double)(boost::posix_time::microsec_clock::local_time() - time0).total_microseconds() / 1000 << '\n';
						__PDEBUG(std::cout << "Done polynomial hash coded multiplying\n");
						// Decode and insert into retval.
						// TODO: add debug info about cms' size here.
						const c_iterator c_it_f = cms.end();
						term_type1 tmp_term;
						for (c_iterator c_it = cms.begin(); c_it != c_it_f; ++c_it) {
							this->decode(c_it->first,c_it->second + 2 * this->m_fast_h.lower(),tmp_term);
							if (!tmp_term.is_canonical(args_tuple)) {
								tmp_term.canonicalise(args_tuple);
							}
							this->m_retval.insert(tmp_term, args_tuple);
						}
						__PDEBUG(std::cout << "Done polynomial hash coded\n");
					}
			};
	};
}

#endif
