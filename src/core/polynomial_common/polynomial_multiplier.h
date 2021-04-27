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


#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_multiplier.h"
#include "../coded_hash_table.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../settings.h" // For debug and cache size.
#include "../stats.h"
#include "../type_traits.h"
#include "../utils.h" // For iota.

#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/thread/thread.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/integral_constant.hpp>

#include <cstddef>
#include <exception>
#include <stdexcept>
#include <utility> // For std::pair.
#include <vector>
#include <algorithm> // For std::max.


namespace piranha
{
// Generic threaded vector multiplication.
template <class Functor>
struct ThreadedBlockedMultiplier
{
	ThreadedBlockedMultiplier(const std::size_t &nominal_block_size, const std::size_t &size1, const std::size_t &size2, const std::size_t &thread_id,
		                        const std::size_t &thread_n, boost::barrier *barrier, std::size_t &cur_idx1_start, bool &breakout, Functor &func,
		                        BlockSequence &idx_vector1, BlockSequence &idx_vector2)
    : m_nominal_block_size(nominal_block_size), m_size1(size1), m_thread_id(thread_id), m_thread_n(thread_n),
	  m_barrier(barrier), m_cur_idx1_start(cur_idx1_start), m_breakout(breakout), m_func(func), m_idx_vector1(idx_vector1), m_idx_vector2(idx_vector2)
	{
#ifdef _DEBUG
		std::cout << "threaded_blocked_multiplier" << std::endl 
			      << "nominal_block_size: " << nominal_block_size << std::endl 
				  << "size1: " << size1 << std::endl
				  << "size2: " << size2 << std::endl
				  << "thread_id: " << thread_id << std::endl
				  << "thread_n:  " << thread_n << std:: endl
				  << "cur_idx1_start: " << cur_idx1_start <<std::endl
				  << "breakout: " << breakout << std::endl
				  <<std::flush;
#endif
		// Sanity checks.
		PIRANHA_ASSERT(thread_n > 0 && thread_id < thread_n && (barrier || thread_n == 1));

		// Numerical limits check. We need an extra block size buffer at the end to make sure we are able to
		// represent all indices and sizes.
		// TODO: we need to take care of extra thread_n blocks at the end, they must be representable. Maybe not, after all.
		if (nominal_block_size > std::numeric_limits<std::size_t>::max() || size1 >= std::numeric_limits<std::size_t>::max() - nominal_block_size ||
			size2 >= std::numeric_limits<std::size_t>::max() - nominal_block_size)
		{
			PIRANHA_THROW(std::overflow_error, "numerical overflow in threaded block multiplication");
		}
	}


	void operator()()
	{
//		std::cout << "polynomial_multiplier::(): 1 : threadid = " << m_thread_id << std::endl << std::flush;

		if (m_thread_id == 0)
		{
//			std::cout << "polynomial_multiplier::(): 2 : threadid = " << m_thread_id << std::endl << std::flush;
			// The first thread is in charge of the initial setup of the indices vectors.
			// TODO: exception handling, in case of both single and multi thread.
			m_idx_vector1.resize(boost::numeric_cast<BlockSequence::size_type>(m_thread_n));
			m_idx_vector2.resize(boost::numeric_cast<BlockSequence::size_type>(m_thread_n));
			m_cur_idx1_start = 0;
		}

		sync();
//		std::cout << "polynomial_multiplier::(): 3 : threadid = " << m_thread_id << std::endl << std::flush;
		BlockSequence orig2(m_idx_vector2.size());
		
		while (m_cur_idx1_start != m_size1)
		{
//			std::cout << "polynomial_multiplier::(): 4 : threadid = " << m_thread_id << std::endl << std::flush;
//			std::cout << "polynomial_multiplier::(): " << m_cur_idx1_start <<std::endl <<std::flush;
			sync();
			
			if (m_thread_id == 0)
			{
//				std::cout << "polynomial_multiplier::(): 5 : threadid = " << m_thread_id << std::endl << std::flush;
				m_func.blocksSetup(m_cur_idx1_start ,m_nominal_block_size, m_idx_vector1, m_idx_vector2);
				m_breakout = false;
			}

			sync();
//			std::cout << "polynomial_multiplier::(): 6 : threadid = " << m_thread_id << std::endl << std::flush;
			const std::size_t i_start = m_idx_vector1[m_thread_id].first, i_end = m_idx_vector1[m_thread_id].second;
			// Remember the original block sequence for the second series.
			std::copy(m_idx_vector2.begin(), m_idx_vector2.end(), orig2.begin());
			// Reset the wrap count.
			std::size_t wrap_count = 0;

			while (true)
			{
//				std::cout << "polynomial_multiplier::(): 7 : threadid = " << m_thread_id << std::endl << std::flush;;
				const std::size_t j_start = m_idx_vector2[m_thread_id].first, j_end = m_idx_vector2[m_thread_id].second;
				// NOTE: here maybe we can put a preemptive check on j start/end so that if the inner
				// cycle is empty we skip this part altogether.
				for (std::size_t i = i_start; i < i_end; ++i)
				{
					for (std::size_t j = j_start; j < j_end; ++j)
					{
						// TODO: truncation.
						m_func(i,j);

					}
				}

				sync();
//				std::cout << "polynomial_multiplier::(): 8 : threadid = " << m_thread_id << std::endl << std::flush;

				if (m_thread_id == 0)
				{
//					std::cout << "polynomial_multiplier::(): 9 : threadid = " << m_thread_id << std::endl << std::flush;
					if (!m_func.block2_advance(m_idx_vector1,m_idx_vector2,m_nominal_block_size,orig2,wrap_count))
					{
//						std::cout << "polynomial_multiplier::(): 9a : break" << std::endl << std::flush;
						m_breakout = true;
					}
				}
//				std::cout << "polynomial_multiplier::(): 10a : threadid = " << m_thread_id << std::endl << std::flush;
				sync();
//				std::cout << "polynomial_multiplier::(): 10 : threadid = " << m_thread_id << std::endl << std::flush;

				if (m_breakout)
				{
					break;
				}
			}
		}
	}


	void sync()
	{
		if (m_barrier)
		{
//			std::cout << "polynomial_multiplier::sync(): threadid = " << m_thread_id << std::endl << std::flush;
			m_barrier->wait();
//			std::cout << "polynomial_multiplier::sync(): finished : threadid = " << m_thread_id << std::endl << std::flush;
		}
	}


	const std::size_t	m_nominal_block_size;
	const std::size_t	m_size1;
	const std::size_t	m_thread_id;
	const std::size_t	m_thread_n;
	boost::barrier	   *m_barrier;
	std::size_t		   &m_cur_idx1_start;
	bool			   &m_breakout;
	Functor			   &m_func;
	BlockSequence	   &m_idx_vector1;
	BlockSequence	   &m_idx_vector2;
};


template <class Series1, class Series2, class ArgsTuple, class GenericTruncator>
struct PolynomialVectorFunctor: public BaseCodedFunctor<Series1, Series2, ArgsTuple, GenericTruncator, PolynomialVectorFunctor<Series1, Series2, ArgsTuple, GenericTruncator> >
{
	typedef typename FinalCf<Series1>::Type cf_type1;
	typedef typename FinalCf<Series2>::Type cf_type2;
	typedef typename Series1::TermType term_type1;
	typedef typename Series2::TermType term_type2;
	typedef BaseCodedFunctor<Series1, Series2, ArgsTuple, GenericTruncator, PolynomialVectorFunctor<Series1, Series2, ArgsTuple, GenericTruncator> > Ancestor;

	PolynomialVectorFunctor(std::vector<cf_type1>            &tc1,       std::vector<cf_type2> &tc2,
		                    std::vector<MaxFastInt>        &ck1,       std::vector<MaxFastInt> &ck2,
		                    std::vector<term_type1 const *>  &t1,        std::vector<term_type2 const *> &t2,
		                    const GenericTruncator           &truncator, cf_type1 *vc_res, const ArgsTuple &argsTuple)
        : Ancestor(tc1, tc2, ck1, ck2, t1, t2, truncator, argsTuple), m_vc_res(vc_res)
	{}


	bool operator()(std::size_t const i, std::size_t const j)
	{
		if (this->m_trunc.skip(&this->m_t1[i], &this->m_t2[j]))
        {
			return false;
		}
		// Calculate index of the result.
		const MaxFastInt res_index = this->m_ck1[i] + this->m_ck2[j];
		m_vc_res[res_index].addmul(this->m_tc1[i], this->m_tc2[j], this->m_argsTuple);

		return true;
	}


	void blocksSetup(std::size_t & currentIndexStart, std::size_t const blockSize, BlockSequence &indexVector1, BlockSequence &indexVector2)
	{
		if (currentIndexStart == 0)
		{
			initial_setup();
		}
		this->baseBlocksSetup(currentIndexStart, blockSize, indexVector1, indexVector2);
	}


	static const MaxFastInt &get_mem_pos(const MaxFastInt &n)
	{
		return n;
	}


	// Return the two intervals in indices in the output structure containing the results of
	// one block-by-block multiplication.
	std::pair<BlockInterval, BlockInterval> blocksToIntervals(BlockType const &b1, BlockType const &b2) const
	{
		PIRANHA_ASSERT(b1.first <= b1.second && b2.first <= b2.second);

		// If at least one of the blocks is empty, then we won't be writing into any interval of indices.
		if (b1.first == b1.second || b2.first == b2.second)
        {
			return std::make_pair(BlockInterval::empty(), BlockInterval::empty());
		}

		// In case of vector coded, we always end up with a single interval in output.
		return std::make_pair(BlockInterval(this->m_ck1[b1.first] + this->m_ck2[b2.first], this->m_ck1[b1.second - 1] + this->m_ck2[b2.second - 1]),
			                  BlockInterval::empty());
	}


	void initial_setup()
	{
		// Build the permutation vectors.
		typedef std::vector<std::size_t>::size_type size_type;
		
        std::vector<std::size_t> perm1(boost::numeric_cast<size_type>(this->m_ck1.size()));
		std::vector<std::size_t> perm2(boost::numeric_cast<size_type>(this->m_ck2.size()));
		iota(perm1.begin(), perm1.end(), std::size_t(0));
		iota(perm2.begin(), perm2.end(), std::size_t(0));

		// Sort the permutation vectors.
		typedef typename Ancestor::template IndirectSorter<PolynomialVectorFunctor> IndirectSorter;
		
        std::sort(perm1.begin(), perm1.end(), IndirectSorter(*this, this->m_ck1));
		std::sort(perm2.begin(), perm2.end(), IndirectSorter(*this, this->m_ck2));

		// Apply the permutations to the other vectors.
		this->apply_permutation(perm1, this->m_tc1);
		this->apply_permutation(perm1, this->m_ck1);
		this->apply_permutation(perm1, this->m_t1);
		this->apply_permutation(perm2, this->m_tc2);
		this->apply_permutation(perm2, this->m_ck2);
		this->apply_permutation(perm2, this->m_t2);
	}


	cf_type1 *m_vc_res;
};


template <class Series1, class Series2, class ArgsTuple, class GenericTruncator>
struct PolynomialHashFunctor: public BaseCodedFunctor<Series1, Series2, ArgsTuple, GenericTruncator, PolynomialHashFunctor<Series1, Series2, ArgsTuple, GenericTruncator> >
{
	typedef typename FinalCf<Series1>::Type cf_type1;
	typedef typename FinalCf<Series2>::Type cf_type2;
	typedef typename Series1::TermType term_type1;
	typedef typename Series2::TermType term_type2;
	typedef std::pair<cf_type1, MaxFastInt> cterm_type;
	typedef BaseCodedFunctor<Series1, Series2, ArgsTuple, GenericTruncator, PolynomialHashFunctor<Series1, Series2, ArgsTuple, GenericTruncator> > Ancestor;
	typedef coded_hash_table<cf_type1, MaxFastInt, CountingAllocator<char> > csht_type;

	PolynomialHashFunctor(cterm_type &cterm, std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
		std::vector<MaxFastInt> &ck1, std::vector<MaxFastInt> &ck2,
		std::vector<const term_type1 *> &t1, std::vector<const term_type2 *> &t2,
		const GenericTruncator &trunc, csht_type *cms, const ArgsTuple &argsTuple)
    :  Ancestor(tc1, tc2, ck1, ck2, t1, t2, trunc, argsTuple), m_cterm(cterm), m_cms(cms)/*,m_overflow_terms()*/
	{}


	bool operator()(std::size_t const i, std::size_t const j)
	{
		typedef typename csht_type::iterator c_iterator;

		if (this->m_trunc.skip(&this->m_t1[i], &this->m_t2[j])) 
		{
			return false;
		}

		m_cterm.second  = this->m_ck1[i];
		m_cterm.second += this->m_ck2[j];
		std::pair<bool, c_iterator> res = m_cms->find(m_cterm.second);
		
        if (res.first)
		{
			res.second->first.addmul(this->m_tc1[i], this->m_tc2[j], this->m_argsTuple);
		} else
		{
			// Assign to the temporary term the old cf (new_key is already assigned).
			m_cterm.first = this->m_tc1[i];
			// Multiply the old term by the second term.
			m_cterm.first.multBy(this->m_tc2[j], this->m_argsTuple);
			m_cms->insert_new(m_cterm, res.second);
		}

		return true;
	}


	MaxFastInt get_mem_pos(const MaxFastInt &n) const
	{
		// Assert on the convertibility to MaxFastInt. We should be sure because of the logic in coded_hash_table,
		// but just in case...
		PIRANHA_ASSERT(boost::numeric_cast<MaxFastInt>(m_cms->get_memory_position(n)) >= 0);
		return m_cms->get_memory_position(n);
	}


	// Return the two intervals in indices in the output structure containing the results of
	// one block-by-block multiplication.
	std::pair<BlockInterval, BlockInterval> blocks_to_intervals(BlockType const &b1, BlockType const &b2) const
	{
		PIRANHA_ASSERT(b1.first <= b1.second && b2.first <= b2.second);

		// If at least one of the blocks is empty, then we won't be writing into any interval of indices.
		if (b1.first == b1.second || b2.first == b2.second)
		{
			return std::make_pair(BlockInterval::empty(), BlockInterval::empty());
		}

		BlockInterval tmp(get_mem_pos(this->m_ck1[b1.first])      + get_mem_pos(this->m_ck2[b2.first]),
			              get_mem_pos(this->m_ck1[b1.second - 1]) + get_mem_pos(this->m_ck2[b2.second - 1]));

		const MaxFastInt ht_size = boost::numeric_cast<MaxFastInt>(m_cms->get_vector_size());

		PIRANHA_ASSERT(ht_size > 0 && tmp.lower() >= 0 && tmp.upper() >= 0);

		if (tmp.lower() >= ht_size || tmp.upper() < ht_size)
		{
			PIRANHA_ASSERT(tmp.upper() / 2 < ht_size);
			return std::make_pair(BlockInterval(tmp.lower() % ht_size,tmp.upper() % ht_size), BlockInterval::empty());
		} else
		{
			return std::make_pair(BlockInterval(tmp.lower(), ht_size - 1), BlockInterval(0, tmp.upper() % ht_size));
		}
	}


	cterm_type		&m_cterm;
	csht_type		*m_cms;
	//std::vector<cterm_type>	&m_overflow_terms;
};


/// Series multiplier specifically tuned for polynomials.
/**
 * This multiplier internally will use coded arithmetics if possible, otherwise it will operate just
 * like piranha::BaseSeriesMultiplier.
 */
struct polynomial_multiplier
{
	template <class Series1, class Series2, class ArgsTuple, class Truncator>
	class get_type: public BaseSeriesMultiplier< Series1, Series2, ArgsTuple, Truncator, get_type<Series1, Series2, ArgsTuple, Truncator> >,
		            public CodedMultiplier<get_type<Series1, Series2, ArgsTuple, Truncator>, Series1, Series2, boost::tuple<boost::true_type> >
	{
			typedef BaseSeriesMultiplier< Series1, Series2, ArgsTuple, Truncator, get_type<Series1, Series2, ArgsTuple, Truncator> >    ancestor;
			typedef CodedMultiplier<get_type<Series1, Series2, ArgsTuple, Truncator>, Series1, Series2, boost::tuple<boost::true_type> > coded_ancestor;
			
			friend class CodedMultiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,boost::tuple<boost::true_type> >;

			typedef typename ancestor::TermType1 term_type1;
			typedef typename ancestor::TermType2 term_type2;
			typedef typename term_type1::CfType  cf_type1;
			typedef typename term_type2::CfType  cf_type2;

		public:

			typedef Series1   series_type1;
			typedef Series2   series_type2;
			typedef ArgsTuple ArgsTupleType;
			typedef typename Truncator::template GetType<Series1,Series2,ArgsTuple> truncator_type;

			get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &argsTuple)
                : ancestor(s1, s2, retval, argsTuple) {}


			template <class GenericTruncator>
			bool performVectorCodedMultiplication(std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
				                                  std::vector<term_type1 const *> &t1, std::vector<term_type2 const *> &t2, const GenericTruncator &trunc)
			{
				std::vector<cf_type1, CountingAllocator<cf_type1> > vc;
				// Try to allocate the space for vector coded multiplication.
				// The +1 is needed because we need the number of possible codes between min and max.
				PIRANHA_ASSERT(boost::numeric::width(this->m_fast_h) + 1 >= 0);
				const std::size_t n_codes = boost::numeric_cast<std::size_t>(boost::numeric::width(this->m_fast_h) + 1);
				try {
					vc.resize(n_codes);
				} catch (const std::bad_alloc &)
				{
					PIRANHA_DEBUG(std::cout << "Not enough physical memory available for vector coded.\n");

					return false;
				} catch (const memory_error &)
				{
					PIRANHA_DEBUG(std::cout << "Memory limit reached for vector coded.\n");

					return false;
				}
				PIRANHA_DEBUG(std::cout << "Going for vector coded polynomial multiplication\n");
				// Define the base pointers for storing the results of multiplication.
				// NOTE: even if here it seems like we are going to write outside allocated memory,
				//       the indices from the analysis of the coded series will prevent out-of-boundaries
				//       reads/writes. The thing works like this: we have ncodes slots allocated, so memory
				//       indices in [0,ncodes - 1]. But since we are doing arithmetics on shifted codes, the
				//       code range is [-h_min,ncodes - 1 - h_min]. So we shift the baseline memory location
				//       so that we can use shifted codes directly as indices.
				const std::size_t size1 = this->terms1.size();
				const std::size_t size2 = this->terms2.size();
 // sizes std::cout << "sizes: " << size1 << ',' << size2 << '\n';

				PIRANHA_ASSERT(size1 && size2);

                const ArgsTupleType &argsTuple = this->argsTuple;
				cf_type1 *vc_res =  &vc[0] - this->m_fast_h.lower();

				// Find out a suitable block size.
				const std::size_t block_size = this->template computeBlockSize<sizeof(cf_type1)>();

				PIRANHA_DEBUG(std::cout << "Block size: " << block_size << '\n');
				PIRANHA_DEBUG(std::cout << "Block size: " << block_size << '\n');
				
                // Perform multiplication.
				typedef PolynomialVectorFunctor<Series1, Series2, ArgsTuple, GenericTruncator> VectorFunctorType;
				VectorFunctorType vm(tc1, tc2, this->m_ckeys1, this->m_ckeys2a, t1, t2, trunc, vc_res, argsTuple);
//				const std::size_t nthread = settings::get_nthread();
				//TODO:GUT corrected below. There are problems with the number of threads in several places. This is one.
				const std::size_t nthread = std::min(settings::get_nthread(), std::min(size1, size2));
	const boost::posix_time::ptime time0 = boost::posix_time::microsec_clock::local_time();
				
                // Variables needed by the multiplier.
				BlockSequence s1;
                BlockSequence s2;
				bool breakout = false;
				std::size_t cur_idx1_start = 0;
				// TODO: probably we need to rethink a bit this, taking into account also the number of threads. Also drop the truncation limitation.
				if (trunc.isEffective() || (this->terms1.size() * this->terms2.size()) <= 400 || nthread == 1) 
                {
					stats::trace_stat("mult_st", std::size_t(0), increment);
					ThreadedBlockedMultiplier<VectorFunctorType> t(block_size, size1, size2, 0, 1, 0, cur_idx1_start, breakout, vm, s1, s2);
					t();
				} else 
                {
					PIRANHA_DEBUG(std::cout << "using " << nthread << " threads\n");
					stats::trace_stat("mult_mt", std::size_t(0), increment);
					boost::thread_group tg;
					boost::barrier b(static_cast<unsigned int>(nthread));
					for (std::size_t i = 0; i < nthread; ++i) 
					{
						tg.create_thread(ThreadedBlockedMultiplier<VectorFunctorType>(block_size, size1, size2, i, nthread, &b, cur_idx1_start, breakout, vm, s1, s2));
					}
					tg.join_all();
				}
		std::cout << "Elapsed time: " << (double)(boost::posix_time::microsec_clock::local_time() - time0).total_microseconds() << " micro seconds" <<std::endl;
				PIRANHA_DEBUG(std::cout << "Done multiplying\n");

				const MaxFastInt i_f = this->m_fast_h.upper();
				// Decode and insert the results into return value.
				term_type1 tmp_term;
				for (MaxFastInt i = this->m_fast_h.lower(); i <= i_f; ++i)
				{
					// Take a shortcut and check for ignorability of the coefficient here.
					// This way we avoid decodification, and all the series term insertion yadda-yadda.
					if (!vc_res[i].isIgnorable(argsTuple))
					{
						this->decode(vc_res[i], i, tmp_term);
						if (!tmp_term.isCanonical(argsTuple))
						{
							tmp_term.canonicalise(argsTuple);
						}
						this->retval.insert(tmp_term, argsTuple);
					}
				}
				PIRANHA_DEBUG(std::cout << "Done polynomial vector coded.\n");
				return true;
			}


			template <class GenericTruncator>
			void performHashCodedMultiplication(std::vector<cf_type1>           &tc1, std::vector<cf_type2> &tc2,
				                                std::vector<const term_type1 *> &t1,  std::vector<const term_type2 *> &t2,
                                                GenericTruncator const &truncator)
			{
				typedef coded_hash_table<cf_type1, MaxFastInt, CountingAllocator<char> > csht;

				typedef typename csht::iterator c_iterator;
				stats::trace_stat("mult_st", std::size_t(0), increment);
				// Let's find a sensible size hint.
				const std::size_t n_codes = boost::numeric_cast<std::size_t>(boost::numeric::width(this->m_fast_h) + 1);
				const std::size_t size_hint = static_cast<std::size_t>(
					std::max<double>(this->m_density1,this->m_density2) * n_codes);
				const std::size_t size1 = this->terms1.size(); 
                const std::size_t size2 = this->terms2.size();

				PIRANHA_ASSERT(size1 && size2);

				const ArgsTupleType &argsTuple = this->argsTuple;
				csht cms(size_hint);
				// Find out a suitable block size.
				const std::size_t block_size = this->template computeBlockSize<sizeof(std::pair<cf_type1, MaxFastInt>)>();

				PIRANHA_DEBUG(std::cout << "Block size: " << block_size << '\n');

                //std::cout << "Block size: " << block_size << '\n';
                const boost::posix_time::ptime time0 = boost::posix_time::microsec_clock::local_time();
				
                std::pair<cf_type1, MaxFastInt> cterm;
				PolynomialHashFunctor<Series1, Series2, ArgsTuple, GenericTruncator> hm(cterm, tc1, tc2, this->m_ckeys1, this->m_ckeys2a, t1, t2, truncator, &cms, argsTuple);

				this->blockedMultiplication(block_size, size1, size2, hm);
                
                //std::cout << "Elapsed time: " << (double)(boost::posix_time::microsec_clock::local_time() - time0).total_microseconds() / 1000 << '\n';

				PIRANHA_DEBUG(std::cout << "Done polynomial hash coded multiplying\n");

                // Decode and insert into retval.
				// TODO: add debug info about cms' size here.
				const c_iterator c_it_f = cms.end();
				term_type1 tmp_term;

				for (c_iterator c_it = cms.begin(); c_it != c_it_f; ++c_it)
				{
					this->decode(c_it->first, c_it->second + 2 * this->m_fast_h.lower(), tmp_term);
					if (!tmp_term.isCanonical(argsTuple)) 
					{
						tmp_term.canonicalise(argsTuple);
					}
					this->retval.insert(tmp_term, argsTuple);
				}

				PIRANHA_DEBUG(std::cout << "Done polynomial hash coded\n");
			}
	};
};
}

#endif
