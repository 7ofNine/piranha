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



#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"
#include "../settings.h"
#include "../stats.h"
#include "../type_traits.h"
#include "../lambdas.h"
#include "coded_multiplier_mp.h"
#include "null_truncator.h"

#include <boost/functional/hash.hpp>
#include <boost/iterator/permutation_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp> // We assert equality between vh tuples below.

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <type_traits>

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

// TODO:
// - use power of 2 coding for hash coded?
// - cache usage: determine optimal size at runtime, e.g., inspecting the size of MP coefficients?
// - numeric cast in coded functor & friends.


namespace piranha
{

typedef std::pair<std::size_t, std::size_t> BlockType;
typedef std::vector<BlockType>              BlockSequence;
typedef boost::numeric::interval< MaxFastInt, boost::numeric::interval_lib::policies< boost::numeric::interval_lib::rounded_math<MaxFastInt>,
	                                          boost::numeric::interval_lib::checking_base<MaxFastInt> > > 
                                            BlockInterval;


template <class Series1, class Series2, class ArgsTuple, class GenericTruncator, class Derived>
struct BaseCodedFunctor
{
	typedef typename FinalCf<Series1>::Type CfType1;
	typedef typename FinalCf<Series2>::Type CfType2;
	typedef typename Series1::TermType TermType1;
	typedef typename Series2::TermType TermType2;


    // 
    // Functorial class providing the binary predicate (comaprison) operation for a std::sort
    // It is indirect by sorting the indices not the vector itself
	template <class Functor>
	class IndirectSorter
	{
        public:

		IndirectSorter(Functor const &functor, std::vector<MaxFastInt> const &vec) : functor(functor), vec(vec) {}

		bool operator()(std::size_t const n1, std::size_t const n2) const
		{
			// TODO: numeric casts here, or maybe one large check at the beginning of the coded multiplier?
			return functor.get_mem_pos(vec[n1]) < functor.get_mem_pos(vec[n2]);
		}

        private:
		Functor const			        &functor;
		std::vector<MaxFastInt> const	&vec;
	};


	BaseCodedFunctor(std::vector<CfType1>           &tc1, std::vector<CfType2>           &tc2,
		             std::vector<MaxFastInt>      &ck1, std::vector<MaxFastInt>      &ck2,
		             std::vector<TermType1 const *> &t1,  std::vector<TermType2 const *> &t2,
		             GenericTruncator const &truncator,    ArgsTuple const &argsTuple)
    : m_tc1(tc1), m_tc2(tc2), m_ck1(ck1), m_ck2(ck2), m_t1(t1), m_t2(t2), m_trunc(truncator), m_argsTuple(argsTuple)
	{}


	void baseBlocksSetup(std::size_t & currentIndex1Start, std::size_t const blockSize, BlockSequence &indexVector1, BlockSequence &indexVector2)
	{
//#ifdef _DEBUG
//		std::cout << "BaseCodedFunctor::baseBlocksSetup" <<std::endl
//			      << "currentIndex1Start: "              << currentIndex1Start  << std::endl
//				  << "blockSize: "                       << blockSize           << std::endl
//				  << "indexVector1 size: "               << indexVector1.size() << std::endl
//				  << "indexVector2 size: "               << indexVector2.size() << std::endl;
//#endif
		PIRANHA_ASSERT(currentIndex1Start < m_tc1.size() && indexVector1.size() == indexVector2.size());

		std::size_t upperBound1 = std::min<std::size_t>(m_tc1.size(), currentIndex1Start + indexVector1.size() * blockSize);
		std::size_t	upperBound2 = std::min<std::size_t>(m_tc2.size(), indexVector2.size() * blockSize);

		do {
			// Now we must check the blocks for the following conditions:
			// 1 - we must not be past the end of the series.
			// 2 - the macroblocks must not result in overlapping areas in the output structure.
			// 3 - the upper bound of each block must be different from the lower bound of next block.
			// ---
			// Determine the sequences, given upper and lower bounds.
			determineSequence(indexVector1, currentIndex1Start, upperBound1);
			determineSequence(indexVector2,                  0, upperBound2);

			// Prepare lower-upper bounds for the next iteration, if any. But we never want to have less than 1 term in the whole sequence.
			PIRANHA_ASSERT(upperBound1 > currentIndex1Start);

			if (indexVector1.back().second - currentIndex1Start >= 2) 
			{
				upperBound1 = currentIndex1Start + (indexVector1.back().second - currentIndex1Start) / 2;
			}

			if (indexVector2.back().second >= 2) 
			{
				upperBound2 = indexVector2.back().second / 2;
			}

		} while (sequencesOverlap(indexVector1, indexVector2));

		// Blocks boundaries check.
		adjustBlockBoundaries(indexVector1, m_ck1);
		adjustBlockBoundaries(indexVector2, m_ck2);

		// Finally, update the currentIndex1Start.
		currentIndex1Start = indexVector1.back().second;

// std::cout << "init\n";
// for (std::size_t i = 0; i < indexVector1.size(); ++i) {
// 	std::cout << indexVector1[i].first << ',' << indexVector1[i].second << '\n';
// }
// for (std::size_t i = 0; i < indexVector2.size(); ++i) {
// 	std::cout << indexVector2[i].first << ',' << indexVector2[i].second << '\n';
// }
// std::cout << "blah\n";
	}


	// Write into s a sequence of blocks ranging from index lower_bound to index upper_bound.
	static void determineSequence(BlockSequence &s, std::size_t const lowerBound, std::size_t const upperBound)
	{
		PIRANHA_ASSERT(upperBound > lowerBound && s.size() > 0);
		const std::size_t numBlocks = boost::numeric_cast<std::size_t>(s.size());
        const std::size_t numTerms  = upperBound - lowerBound;
		const std::size_t blockSize = (numBlocks > numTerms) ? 1 : numTerms / numBlocks;
		for (std::size_t i = 0; i < numBlocks - 1; ++i) 
		{
			s[i].first  = std::min<std::size_t>(upperBound, lowerBound +  i      * blockSize);
			s[i].second = std::min<std::size_t>(upperBound, lowerBound + (i + 1) * blockSize);
		}

		// Handle the last block separately, as it might be non homogeneous.
		s[numBlocks - 1].first = std::min<std::size_t>(upperBound, lowerBound + (numBlocks - 1) * blockSize);
		s[numBlocks - 1].second = upperBound;
	}


	void adjustBlockBoundaries(BlockSequence &s, std::vector<MaxFastInt> const &ck) const
	{
//#ifdef _DEBUG
//		std::cout << "BaseCodedFunctor::adjust_block_boundaries: ck size = " << ck.size() << std::endl;
//#endif
		PIRANHA_ASSERT(s.size() > 0);

		for (BlockSequence::size_type i = 0; i < s.size() - 1; ++i) 
		{
			// Shrink non-empty blocks whose upper value's memory position is shared with the next block.
//#ifdef _DEBUG
//			std::cout << "BaseCodedFunctor::adjust_block_boundaries: blockSequence index = " << i << std::endl << std::flush;
//			std::cout << "s[i] = " << s[i].first << "/" << s[i].second << std::endl << std::flush;
//			std::cout << " upper/lower = " << s[i].second - 1 <<"/" << s[i + 1].first << std::endl << std::flush;  
//#endif
			while (s[i].first != s[i].second && (s[i+1].first < ck.size()) && derived_const_cast->get_mem_pos(ck[s[i].second - 1]) ==
				derived_const_cast->get_mem_pos(ck[s[i + 1].first]))
			{
				--s[i].second;
				--s[i + 1].first;
//#ifdef _DEBUG
//			std::cout << "s[i] = " << s[i].first << "/" << s[i].second << std::endl << std::flush;
//			std::cout << " upper/lower = " << s[i].second - 1 <<"/" << s[i + 1].first << std::endl << std::flush;  
//#endif
			}
		}
	}


	bool block2_advance(const BlockSequence &indexVector1, BlockSequence &indexVector2, const std::size_t &BlockSize, const BlockSequence &orig2, std::size_t &wrapCount) const
	{
		PIRANHA_DEBUG(std::cout << "coded_multiplier::block2_advance()" << std::endl << std::flush);
		PIRANHA_ASSERT(indexVector1.size() == indexVector2.size() && indexVector1.size() > 0);

		if (wrapCount)
		{
			PIRANHA_ASSERT(wrapCount < indexVector2.size());
			if (wrapCount == indexVector2.size() - 1)
			{
				// This means we are at the end.
				return false;
			}
			// Shift down the blocks.
			std::copy(indexVector2.begin() + 1, indexVector2.end(), indexVector2.begin());
			// Get the new block from the originals.
			indexVector2.back() = orig2[wrapCount];
			// Increase the wrap count.
			++wrapCount;
		} else
		{
			// Shift down the blocks.
			std::copy(indexVector2.begin() + 1, indexVector2.end(), indexVector2.begin());
			// Set the new starting point for the last block.
			indexVector2.back().first = indexVector2.back().second;
			// Add the block size or stop at the end of the series, if necessary.
			indexVector2.back().second = std::min<std::size_t>(m_tc2.size(), indexVector2.back().first + BlockSize);
			// Now check if we are at the end of the first phase.
			if (indexVector2.front() == BlockType(m_tc2.size(), m_tc2.size()))
			{
				if (indexVector2.size() > 1)
				{
					// If multi-threaded, insert at the end the first original block.
					indexVector2.back() = orig2.front();
					// Start the wrap count.
					wrapCount = 1;
				} else
				{
					// In single-threaded, this means we have finished.
					return false;
				}
			} else if (indexVector2.size() > 1)
			{
				// If we are not at the end of the first phase and we are multithreaded, we need to make sure the newly-added
				// block does not overlap.
				// NOTE: maybe this function can be replaced by direct check that the last block of first series
				// by the newly added block in second series do not overlap with the remaining macroblock 1 by remaining
				// macro block 2.
				while (sequencesOverlap(indexVector1, indexVector2))
				{
					PIRANHA_ASSERT(indexVector2.back().second >= indexVector2.back().first);
					indexVector2.back().second = indexVector2.back().first + (indexVector2.back().second - indexVector2.back().first) / 2;
				}
			}
		}

		// Make sure we have no overlaps.
		PIRANHA_ASSERT(!sequencesOverlap(indexVector1, indexVector2));
// std::cout << "after advance\n";
// for (std::size_t i = 0; i < indexVector1.size(); ++i) {
// 	std::cout << indexVector1[i].first << ',' << indexVector1[i].second << '\n';
// }
// for (std::size_t i = 0; i < indexVector2.size(); ++i) {
// 	std::cout << indexVector2[i].first << ',' << indexVector2[i].second << '\n';
// }
// std::cout << "blappo\n";

		return true;
	}


	static bool intervalSorter(BlockInterval const &i1,  BlockInterval const &i2)
	{
		return i1.lower() < i2.lower();
	}


	bool sequencesOverlap(BlockSequence const &s1,  BlockSequence const &s2) const
	{
		PIRANHA_ASSERT(s1.size() == s2.size() && s1.size() > 0);

		typedef std::vector<BlockInterval>::size_type size_type;

		std::vector<BlockInterval> vi;
		for (size_type i = 0; i < s1.size(); ++i) 
		{
			std::pair<BlockInterval, BlockInterval> tmp(derived_const_cast->blocksToIntervals(s1[i], s2[i]));
			if (!boost::numeric::empty(tmp.first)) 
			{
				vi.push_back(tmp.first);
			}

			if (!boost::numeric::empty(tmp.second)) 
			{
				vi.push_back(tmp.second);
			}
		}

		if (!vi.size()) 
		{
			return false;
		}
		// Sort according to lower bound of the interval.
		std::sort(vi.begin(), vi.end(), intervalSorter);

		PIRANHA_ASSERT(vi.size() > 0);

		// Check that all intervals are disjoint.
		for (size_type i = 0; i < vi.size() - 1; ++i) 
		{
			if (vi[i].upper() >= vi[i + 1].lower()) 
			{
				return true;
			}
		}
		return false;
	}


	// TODO: rewrite with iterators for genericity? Or maybe provide alternative version.
	template <class T>
	static void apply_permutation(const std::vector<std::size_t> &perm, std::vector<T> &v)
	{
		typedef boost::permutation_iterator<typename std::vector<T>::iterator, std::vector<std::size_t>::const_iterator> perm_iterator;
		std::vector<T> other(v.size());
		std::copy(perm_iterator(v.begin(), perm.begin()), perm_iterator(v.end(), perm.end()), other.begin());
		other.swap(v);
	}


	std::vector<CfType1>		    &m_tc1;
	std::vector<CfType2>		    &m_tc2;
	std::vector<MaxFastInt>	    &m_ck1;
	std::vector<MaxFastInt>	    &m_ck2;
	std::vector<TermType1 const *>  &m_t1;
	std::vector<TermType2 const *>  &m_t2;
	const GenericTruncator		    &m_trunc;
	const ArgsTuple			        &m_argsTuple;
};


	/// Toolbox for coded series multiplication.
	/**
	 * Intended to be inherited together with piranha::BaseSeriesMultiplier. It adds common methods for
	 * series multiplication through Kronecker codification. Requirements:
	 * - input series must have at least one term;
	 * - input series must have at least one argument.
	 */
	template <class Derived, class Series1, class Series2, class OpTuple>
	class CodedMultiplier
	{
			// Some static checks.
        static_assert(Series1::echelonLevel == Series2::echelonLevel, "");
        static_assert(boost::tuples::length<OpTuple>::value == Series1::echelonLevel + 1, "");
			// Main typedefs, for internal use.
			// min/max type for input series.
			typedef typename cm_tuple<Series1>::type_minmax minmax_type;
			// multiprecision min/max type.
			typedef typename cm_tuple<Series1>::type_mp_minmax mp_minmax_type;
			// fast min/max type.
			typedef typename cm_tuple<Series1>::type_max_fast_int_minmax fast_minmax_type;
			// Value handler tuples.
			typedef typename cm_tuple<Series1>::type_value_handler value_handler_type;
			// Coding tuples.
			typedef typename cm_tuple<Series1>::type_mp_coding_tuple mp_coding_tuple_type;
			typedef typename cm_tuple<Series1>::type_coding_tuple fast_coding_tuple_type;
			// Decoding tuple.
			typedef typename cm_tuple<Series1>::type_decoding_tuple decoding_tuple_type;
			// These static checks makes sure that the two series have compatible types in the echelon
			// hierarchy, apart from the numerical coefficients.
            static_assert((std::is_same_v<minmax_type, typename cm_tuple<Series2>::type_minmax>), "");
            static_assert((std::is_same_v<value_handler_type, typename cm_tuple<Series2>::type_value_handler>), "");
			
			// Generalised reverse lexicographic comparison.
			class key_revlex_comparison
			{
				public:
			
					template <class Term>
					bool operator()(const Term *t1, const Term *t2) const
					{
						return key_revlex_comparison_impl<Term>::run(t1, t2);
					}
			};


			enum MultiplicationType
			{
				MULTIPLICATION_PLAIN   = 0,
				MULTIPLICATION_VECTOR = 1,
				MULTIPLICATION_HASH   = 2
			};


			static void trace_mult_type(MultiplicationType const type)
			{
				std::string name = "";
				switch (type)
				{
					case MULTIPLICATION_PLAIN:
						name = "multiplication_plain";
						break;
					case MULTIPLICATION_VECTOR:
						name = "multiplication_vector";
						break;
					case MULTIPLICATION_HASH:
						name = "multiplication_hash";
				}
				stats::trace_stat(name, std::size_t(0), increment);
			}

		public:
			/// Default constructor.
			/**
			 * Initialises viability flag to false and sets up data members for future use. Densities are initialised
			 * to zero and m_mp_gr's and m_fast_gr's sizes are initialised according to the arguments tuple stored in
			 * piranha::BaseSeriesMultiplier.
			 */
			CodedMultiplier()
            : m_gr_is_viable(false), m_mp_h(mp_integer(0)), m_fast_h(0), m_density1(0.), m_density2(0.)
			{
				// NOTE: beware the order of inheritance here, make sure to init BaseSeriesMultiplier before,
				//       otherwise m_argsTuple will be uninitialised here.
				// Initialise the member tuples.
				cm_init_vector_tuple<Series1>(m_mp_gr,   derived_const_cast->argsTuple);
				cm_init_vector_tuple<Series1>(m_fast_gr, derived_const_cast->argsTuple);
				cm_init_vector_tuple<Series1>(m_mp_ct,   derived_const_cast->argsTuple);
				cm_init_vector_tuple<Series1>(m_fast_ct, derived_const_cast->argsTuple);
			}


			/// Perform multiplication and place the result into m_retval.
			void performMultiplication()
			{
				const settings::MultiplicationAlgorithm algo = settings::getMultiplicationAlgorithm();
				std::vector<typename Series1::TermType> f_terms1;
				std::vector<typename Series2::TermType> f_terms2;
				// If echelon level is more than zero we need to flatten out the series.
				if (Series1::echelonLevel)
				{
					f_terms1 = derived_cast->series1.flattenTerms(derived_cast->argsTuple);
					f_terms2 = derived_cast->series2.flattenTerms(derived_cast->argsTuple);
					derived_cast->cacheTermsPointers(f_terms1, f_terms2);
				} else
				{
					// Cache term pointers.
					derived_cast->cacheTermsPointers(derived_cast->series1, derived_cast->series2);
				}
				
				//std::cout << "coded_multiplier::performMultiplication  : 0" << std::endl << std::flush;
                // NOTE: hard coded value of 1000.
				if ((algo == settings::MultiplicationAlgorithm::AUTOMATIC && double(derived_cast->terms1.size()) * double(derived_cast->terms2.size()) < 1000)
					|| algo == settings::MultiplicationAlgorithm::PLAIN)
				{
					derived_cast->performPlainMultiplication();
					trace_mult_type(MULTIPLICATION_PLAIN);
					return;
				}

				// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
				// coefficients and coded keys.
				const typename Derived::truncator_type trunc(derived_cast->terms1, derived_cast->terms2, derived_cast->argsTuple);
				determineViability();
				if (!m_gr_is_viable)
				{
					if (algo == settings::MultiplicationAlgorithm::VECTOR_CODED || algo == settings::MultiplicationAlgorithm::HASH_CODED)
                    {
						PIRANHA_THROW(value_error, "coded multiplication requested, but coded representation is not feasible");
					}
					derived_cast->performPlainMultiplication();
					trace_mult_type(MULTIPLICATION_PLAIN);
					return;
				}

				if (trunc.isEffective())
				{
					derived_cast->llPerformMultiplication(trunc);

				} else
				{
					// Sort input series for better cache usage and multi-threaded implementation.
					//std::cout << "coded_multiplier::performMultiplication  : 1" << std::endl << std::flush;
					std::sort(derived_cast->terms1.begin(), derived_cast->terms1.end(), key_revlex_comparison());
					//std::cout << "coded_multiplier::performMultiplication  : 2" << std::endl << std::flush;
					std::sort(derived_cast->terms2.begin(), derived_cast->terms2.end(), key_revlex_comparison());
					//std::cout << "coded_multiplier::performMultiplication  : 3" << std::endl << std::flush;
					derived_cast->llPerformMultiplication(NullTruncator::template GetType<Series1, Series2, typename Derived::ArgsTupleType>(
						                                                                     derived_cast->terms1, derived_cast->terms2, derived_cast->argsTuple));
					//std::cout << "coded_multiplier::performMultiplication  : 4" << std::endl << std::flush;
				}
			}


			template <class GenericTruncator>
			void llPerformMultiplication(const GenericTruncator &trunc)
			{
				const settings::MultiplicationAlgorithm algo = settings::getMultiplicationAlgorithm();
				
                typedef typename FinalCf<Series1>::Type CfType1;
				typedef typename FinalCf<Series2>::Type CfType2;
				// Code terms.
				// NOTE: it is important to code here since at this point we already have sorted input series,
				//       if necessary.
				code_terms();
				// Cache the coefficients.
				std::vector<CfType1> cf1Cache;
				std::vector<CfType2> cf2Cache;
				std::insert_iterator<std::vector<CfType1> > itCache1(cf1Cache, cf1Cache.begin());
				std::insert_iterator<std::vector<CfType2> > itCache2(cf2Cache, cf2Cache.begin());
				std::transform(derived_cast->terms1.begin(), derived_cast->terms1.end(), itCache1, final_cf_getter<Series1>());
				std::transform(derived_cast->terms2.begin(), derived_cast->terms2.end(), itCache2, final_cf_getter<Series2>());
				
                bool vec_res;
				
                if ((algo == settings::MultiplicationAlgorithm::AUTOMATIC && is_sparse()) || algo == settings::MultiplicationAlgorithm::HASH_CODED) 
				{
					vec_res = false;
				} else 
				{
					//std::cout << "coded_multiplier::llPerformMultiplication  1:  " << std::endl << std::flush;
					vec_res = derived_cast->performVectorCodedMultiplication(cf1Cache, cf2Cache, derived_cast->terms1, derived_cast->terms2, trunc);
					//std::cout << "coded_multiplier::llPerformMultiplication  2:  " << std::endl << std::flush;
				}

				if (!vec_res) 
				{
					if (algo == settings::MultiplicationAlgorithm::VECTOR_CODED) 
					{
						PIRANHA_THROW(value_error, "vector coded multiplication requested, but vector coded representation is infeasible");
					}

					shiftCodes();

					derived_cast->performHashCodedMultiplication(cf1Cache, cf2Cache, derived_cast->terms1, derived_cast->terms2, trunc);
					
                    trace_mult_type(MULTIPLICATION_HASH);
				} else 
				{
					//std::cout << "coded_multiplier::llPerformMultiplication  3:  " << std::endl << std::flush;
					trace_mult_type(MULTIPLICATION_VECTOR);
					//std::cout << "coded_multiplier  4:  " << std::endl << std::flush;
				}
				//std::cout << "coded_multiplier::llPerformMultiplication  5:  " << std::endl << std::flush;
			}


			/// Determine whether the global coded representation is viable or not.
			/**
			 * The m_gr_is_viable flag will be set accordingly after this method is called.
			 */
			void determineViability()
			{
				// Make sure that the series have at least one term.
				PIRANHA_ASSERT(derived_const_cast->terms1.size() > 0 && derived_const_cast->terms2.size() > 0);

				// Declare and init the min/max types for the two series.
				minmax_type t1;
                minmax_type t2;
				cm_init_vector_tuple<Series1>(t1, derived_const_cast->argsTuple);
				cm_init_vector_tuple<Series2>(t2, derived_const_cast->argsTuple);

				// Init and test the first series' tuple.
				typedef typename std::vector<typename Series1::TermType const *>::size_type size_type1;
				const size_type1 size1 = derived_const_cast->terms1.size();
				// Value handler tuples.
				value_handler_type vh1, vh2;
				cm_minmax<minmax_type>::run_init(*derived_const_cast->terms1[0], t1, vh1);
				for (size_type1 i = 1; i < size1; ++i) 
				{
					cm_minmax<minmax_type>::run_test(*derived_const_cast->terms1[i], t1, vh1);
				}

				// Init and test the second series' tuple.
				typedef typename std::vector<typename Series2::TermType const *>::size_type size_type2;
				const size_type2 size2 = derived_const_cast->terms2.size();
				cm_minmax<minmax_type>::run_init(*derived_const_cast->terms2[0], t2, vh2);

				for (size_type2 i = 1; i < size2; ++i) 
				{
					cm_minmax<minmax_type>::run_test(*derived_const_cast->terms2[i], t2, vh2);
				}

				// Now compute the global representation in multiprecision.
				cm_global_minmax<OpTuple>::run(t1, vh1, t2, vh2, m_mp_gr);
				// Assign the global value handler tuple.
				PIRANHA_ASSERT(vh1 == vh2);
				m_vh = vh1;
				// Compute the multiprecision coding tuple.
				compute_mp_coding_tuple(m_mp_ct, m_mp_gr);
				// Compute multiprecision codes range.
				tuple_vector_dot(m_mp_gr, m_mp_ct, m_mp_h);
				// To test whether a representation is viable or not, we need to test for the following things:
				// - m_mp_h must be in the MaxFastInt range;
				// - m_mp_h's width must be within halft MaxFastInt's range (needed for 2*chi shifting).
				// Use lexical cast for max interoperability between numerical types.
				// NOTE: here probably we can reduce greatly the number of memory allocations...
//                auto f = [](const std::size_t x) -> std::size_t { return x + 1; }; // increment functor

				if (boost::numeric::subset(m_mp_h,boost::numeric::interval<mp_integer>(
					boost::lexical_cast<mp_integer>(std::numeric_limits<MaxFastInt>::min()),
					boost::lexical_cast<mp_integer>(std::numeric_limits<MaxFastInt>::max()))) &&
					boost::numeric::width(m_mp_h) <=
					boost::lexical_cast<mp_integer>(std::numeric_limits<MaxFastInt>::max()) / 2)
				{
					// Mark representation as viable.
					m_gr_is_viable = true;
					// Log viability.
                    stats::trace_stat("mult_coded_feasible", std::size_t(0), increment);
				} else 
				{
					stats::trace_stat("mult_coded_unfeasible", std::size_t(0), increment);
				}
			}


			/// Code terms.
			void code_terms()
			{
				PIRANHA_ASSERT(m_gr_is_viable);

				// Downcast multiprecision to fast representation.
				cm_mp_tuple_downcast(m_mp_gr, m_fast_gr);
				cm_mp_tuple_downcast(m_mp_ct, m_fast_ct);
				m_fast_h.assign(boost::lexical_cast<MaxFastInt>(m_mp_h.lower()), boost::lexical_cast<MaxFastInt>(m_mp_h.upper()));
				// Build decoding tuple.
				cm_build_decoding_tuple(m_dt, m_fast_gr);
				// Establish if subtraction is requested or not.
				static const bool sub_requested = op_has_sub<OpTuple>::value;
				// Resize codes vectors.
				typedef std::vector<MaxFastInt>::size_type size_type;
				const size_type csize1 = boost::numeric_cast<size_type>(derived_const_cast->terms1.size());
				const size_type csize2 = boost::numeric_cast<size_type>(derived_const_cast->terms2.size());
				m_ckeys1.resize(csize1);
				m_ckeys2a.resize(csize2);
				if (sub_requested) 
				{
					m_ckeys2b.resize(csize2);
				}
				// Now fill in the codes.
				MaxFastInt code_a = 0, code_b = 0;
				for (size_type i = 0; i < csize1; ++i) 
				{
					cm_code<OpTuple>(m_fast_ct, *derived_const_cast->terms1[i], m_vh,code_a, code_b);
					m_ckeys1[i] = code_a;
				}

				for (size_type i = 0; i < csize2; ++i) 
				{
					cm_code<OpTuple>(m_fast_ct, *derived_const_cast->terms2[i], m_vh,code_a, code_b);
					m_ckeys2a[i] = code_a;
					if (sub_requested) 
					{
						m_ckeys2b[i] = code_b;
					}
				}
				// Compute densities.
				const MaxFastInt w = boost::numeric::width(m_fast_h) + 1;
				m_density1 = static_cast<double>(csize1) / w;
				m_density2 = static_cast<double>(csize2) / w;
			}


			/// Decode.
			/**
			 * Decode given code into return value term, using final_cf as the coefficient at the end of the echelon recursion.
			 */
			template <class FinalCf>
			void decode(const FinalCf &final_cf, const MaxFastInt &code, typename Series1::TermType &term) const
			{
				cm_decode(final_cf, m_dt, m_fast_gr, term, m_vh, code, m_fast_h.lower(), derived_const_cast->argsTuple);
			}


			/// Determine whether coded representation is sparse.
			/**
			 * Must be called only if representation is viable, otherwise runtime assertion will fail. Density is compared
			 * against value hard-coded internally.
			 */
			bool is_sparse() const
			{
				// Magic value established empirically. Possibly subject to tuning in the future.
				static const double limit = 1E-4;
				// We don't want this to be called if we haven't established the suitability
				// of the coded representation first.
				PIRANHA_ASSERT(m_gr_is_viable);
				const double max_density = std::max<double>(m_density1,m_density2);
				return (max_density < limit);
			}


			/// Shift codes.
			/**
			 * Move all the the codes so that the minimum code of the representation is 0.
			 */
			void shiftCodes()
			{
				PIRANHA_ASSERT(m_gr_is_viable);

				typedef std::vector<MaxFastInt>::size_type size_type;
				
                const size_type size1  = m_ckeys1.size();
                const size_type size2a = m_ckeys2a.size();
                const size_type size2b = m_ckeys2b.size();
				const MaxFastInt chi = m_fast_h.lower();

				for (size_type i = 0; i < size1; ++i) 
				{
					m_ckeys1[i] -= chi;
				
                	PIRANHA_ASSERT(m_ckeys1[i] >= 0);
				}

				for (size_type i = 0; i < size2a; ++i) 
				{
					m_ckeys2a[i] -= chi;

					PIRANHA_ASSERT(m_ckeys2a[i] >= 0);
				}

				for (size_type i = 0; i < size2b; ++i) 
				{
					m_ckeys2b[i] -= chi;
				
                	PIRANHA_ASSERT(m_ckeys2b[i] >= 0);
				}
			}

		protected:
			/// Is global coded representation viable?
			bool					m_gr_is_viable;
			/// Multiprecision min/max values for the global representation.
			mp_minmax_type				m_mp_gr;
			/// Fast min/max values for the global representation.
			fast_minmax_type			m_fast_gr;
			/// Global value handler tuple.
			value_handler_type			m_vh;
			/// Multiprecision coding tuple.
			mp_coding_tuple_type			m_mp_ct;
			/// Fast coding tuple.
			fast_coding_tuple_type			m_fast_ct;
			/// Decoding tuple.
			decoding_tuple_type			m_dt;
			/// Multiprecision codes range.
			boost::numeric::interval<mp_integer>	m_mp_h;
			/// Fast codes range.
			boost::numeric::interval<MaxFastInt>	m_fast_h;
			/// Codes for the first series.
			std::vector<MaxFastInt>		m_ckeys1;
			/// Codes for the second series, plus.
			std::vector<MaxFastInt>		m_ckeys2a;
			/// Codes for the second series, minus.
			std::vector<MaxFastInt>		m_ckeys2b;
			/// Density of the first series.
			double					m_density1;
			/// Density of the second series.
			double					m_density2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
