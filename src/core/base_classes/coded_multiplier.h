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

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/integer_traits.hpp>
#include <boost/iterator/permutation_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp> // We assert equality between vh tuples below.
#include <boost/type_traits/is_same.hpp>
#include <cstddef>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"
#include "../settings.h"
#include "../stats.h"
#include "../type_traits.h"
#include "coded_multiplier_mp.h"
#include "null_truncator.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

// TODO:
// - use power of 2 coding for hash coded?
// - cache usage: determine optimal size at runtime, e.g., inspecting the size of MP coefficients?
// - numeric cast in coded functor & friends.

namespace piranha
{

typedef std::pair<std::size_t,std::size_t> block_type;
typedef std::vector<block_type> block_sequence;
typedef boost::numeric::interval<max_fast_int,boost::numeric::interval_lib::policies<
	boost::numeric::interval_lib::rounded_math<max_fast_int>,
	boost::numeric::interval_lib::checking_base<max_fast_int>
> > block_interval;

template <class Series1, class Series2, class ArgsTuple, class GenericTruncator, class Derived>
struct base_coded_functor
{
	typedef typename final_cf<Series1>::type cf_type1;
	typedef typename final_cf<Series2>::type cf_type2;
	typedef typename Series1::term_type term_type1;
	typedef typename Series2::term_type term_type2;
	template <class Functor>
	struct indirect_sorter
	{
		indirect_sorter(const Functor &func, const std::vector<max_fast_int> &v):m_func(func),m_v(v) {}
		bool operator()(const std::size_t &n1, const std::size_t &n2) const
		{
			// TODO numeric casts here, or maybe one large check at the beginning of the coded multiplier?
			return m_func.get_mem_pos(m_v[n1]) < m_func.get_mem_pos(m_v[n2]);
		}
		const Functor			&m_func;
		const std::vector<max_fast_int>	&m_v;
	};
	base_coded_functor(std::vector<cf_type1> &tc1, std::vector<cf_type2> &tc2,
		std::vector<max_fast_int> &ck1, std::vector<max_fast_int> &ck2,
		std::vector<term_type1 const *> &t1, std::vector<term_type2 const *> &t2,
		const GenericTruncator &trunc, const ArgsTuple &args_tuple):
		m_tc1(tc1),m_tc2(tc2),m_ck1(ck1),m_ck2(ck2),m_t1(t1),m_t2(t2),m_trunc(trunc),m_args_tuple(args_tuple)
	{}
	void base_blocks_setup(std::size_t &cur_idx1_start, const std::size_t &block_size,
		block_sequence &idx_vector1, block_sequence &idx_vector2)
	{
		piranha_assert(cur_idx1_start < m_tc1.size() && idx_vector1.size() == idx_vector2.size());
		std::size_t upper_bound1 = std::min<std::size_t>(m_tc1.size(),cur_idx1_start + idx_vector1.size() * block_size),
			upper_bound2 = std::min<std::size_t>(m_tc2.size(),idx_vector2.size() * block_size);
		do {
			// Now we must check the blocks for the following conditions:
			// 1 - we must not be past the end of the series.
			// 2 - the macroblocks must not result in overlapping areas in the output structure.
			// 3 - the upper bound of each block must be different from the lower bound of next block.
			// ---
			// Determine the sequences, given upper and lower bounds.
			determine_sequence(idx_vector1,cur_idx1_start,upper_bound1);
			determine_sequence(idx_vector2,0,upper_bound2);
			// Prepare lower-upper bounds for the next iteration, if any. But we never want to have less than 1 term in the whole sequence.
			piranha_assert(upper_bound1 > cur_idx1_start);
			if (idx_vector1.back().second - cur_idx1_start >= 2) {
				upper_bound1 = cur_idx1_start + (idx_vector1.back().second - cur_idx1_start) / 2;
			}
			if (idx_vector2.back().second >= 2) {
				upper_bound2 = idx_vector2.back().second / 2;
			}
		} while (sequences_overlap(idx_vector1,idx_vector2));
		// 3 - Blocks boundaries check.
		adjust_block_boundaries(idx_vector1,idx_vector2);
		// Finally, update the cur_idx1.
		cur_idx1_start = idx_vector1.back().second;
// std::cout << "init\n";
// for (std::size_t i = 0; i < idx_vector1.size(); ++i) {
// 	std::cout << idx_vector1[i].first << ',' << idx_vector1[i].second << '\n';
// }
// for (std::size_t i = 0; i < idx_vector2.size(); ++i) {
// 	std::cout << idx_vector2[i].first << ',' << idx_vector2[i].second << '\n';
// }
// std::cout << "blah\n";
	}
	// Write into s a sequence of blocks ranging from index lower_bound to at most index upper_bound.
	static void determine_sequence(block_sequence &s, const std::size_t &lower_bound, const std::size_t &upper_bound)
	{
		piranha_assert(upper_bound > lower_bound && s.size() > 0);
		const std::size_t n_blocks = boost::numeric_cast<std::size_t>(s.size()), n_terms = upper_bound - lower_bound;
		const std::size_t block_size = (n_blocks > n_terms) ? 1 : n_terms / n_blocks;
		std::size_t i = 0;
		for (; i < n_blocks - 1; ++i) {
			s[i].first = std::min<std::size_t>(upper_bound,lower_bound + i * block_size);
			s[i].second = std::min<std::size_t>(upper_bound,lower_bound + (i + 1) * block_size);
		}
		// Handle the last block separately, as it might be non homogeneous.
		s[i].first = std::min<std::size_t>(upper_bound,lower_bound + i * block_size);
		s[i].second = upper_bound;
	}
	void adjust_block_boundaries(block_sequence &, block_sequence &)
	{
	}
	bool block2_advance(const block_sequence &idx_vector1, block_sequence &idx_vector2,
		const std::size_t &block_size, const block_sequence &orig2, std::size_t &wrap_count) const
	{
		piranha_assert(idx_vector1.size() == idx_vector2.size() && idx_vector1.size() > 0);
		if (wrap_count) {
			piranha_assert(wrap_count < idx_vector2.size());
			if (wrap_count == idx_vector2.size() - 1) {
				// This means we are at the end.
				return false;
			}
			// Shift down the blocks.
			std::copy(idx_vector2.begin() + 1,idx_vector2.end(),idx_vector2.begin());
			// Get the new block from the originals.
			idx_vector2.back() = orig2[wrap_count];
			// Increase the wrap count.
			++wrap_count;
		} else {
			// Shift down the blocks.
			std::copy(idx_vector2.begin() + 1,idx_vector2.end(),idx_vector2.begin());
			// Set the new starting point for the last block.
			idx_vector2.back().first = idx_vector2.back().second;
			// Add the block size or stop at the end of the series, if necessary.
			idx_vector2.back().second = std::min<std::size_t>(m_tc2.size(),idx_vector2.back().first + block_size);
			// Now check if we are at the end of the first phase.
			if (idx_vector2.front() == block_type(m_tc2.size(),m_tc2.size())) {
				if (idx_vector2.size() > 1) {
					// If multi-threaded, insert at the end the first original block.
					idx_vector2.back() = orig2.front();
					// Start the wrap count.
					wrap_count = 1;
				} else {
					// In single-threaded, this means we have finished.
					return false;
				}
			} else if (idx_vector2.size() > 1) {
				// If we are not at the end of the first phase and we are multithreaded, we need to make sure the newly-added
				// block does not overlap.
				// NOTE: maybe this function can be replaced by direct check that the last block of first series
				// by the newly added block in second series do not overlap with the remaining macroblock 1 by remaining
				// macro block 2.
				while (sequences_overlap(idx_vector1,idx_vector2)) {
					piranha_assert(idx_vector2.back().second >= idx_vector2.back().first);
					idx_vector2.back().second = idx_vector2.back().first +
						(idx_vector2.back().second - idx_vector2.back().first) / 2;
				}
			}
		}
		// Make sure we have no overlaps.
		piranha_assert(!sequences_overlap(idx_vector1,idx_vector2));
// std::cout << "after advance\n";
// for (std::size_t i = 0; i < idx_vector1.size(); ++i) {
// 	std::cout << idx_vector1[i].first << ',' << idx_vector1[i].second << '\n';
// }
// for (std::size_t i = 0; i < idx_vector2.size(); ++i) {
// 	std::cout << idx_vector2[i].first << ',' << idx_vector2[i].second << '\n';
// }
// std::cout << "blappo\n";
		return true;
	}
	static bool interval_sorter(const block_interval &i1, const block_interval &i2)
	{
		return i1.lower() < i2.lower();
	}
	bool sequences_overlap(const block_sequence &s1, const block_sequence &s2) const
	{
		piranha_assert(s1.size() == s2.size() && s1.size() > 0);
		typedef std::vector<block_interval>::size_type size_type;
		std::vector<block_interval> vi;
		for (size_type i = 0; i < s1.size(); ++i) {
			std::pair<block_interval,block_interval> tmp(derived_const_cast->blocks_to_intervals(s1[i],s2[i]));
			if (!boost::numeric::empty(tmp.first)) {
				vi.push_back(tmp.first);
			}
			if (!boost::numeric::empty(tmp.second)) {
				vi.push_back(tmp.second);
			}
		}
		if (!vi.size()) {
			return false;
		}
		// Sort according to lower bound of the interval.
		std::sort(vi.begin(),vi.end(),interval_sorter);
		piranha_assert(vi.size() > 0);
		// Check that all intervals are disjoint.
		for (std::vector<block_interval>::size_type i = 0; i < vi.size() - 1; ++i) {
			if (vi[i].upper() >= vi[i + 1].lower()) {
				return true;
			}
		}
		return false;
	}
	// TODO: rewrite with iterators for genericity? Or maybe provide alternative version.
	template <class T>
	static void apply_permutation(const std::vector<std::size_t> &perm, std::vector<T> &v)
	{
		typedef boost::permutation_iterator<typename std::vector<T>::iterator,std::vector<std::size_t>::const_iterator> perm_iterator;
		std::vector<T> other(v.size());
		std::copy(perm_iterator(v.begin(),perm.begin()),perm_iterator(v.end(),perm.end()),other.begin());
		other.swap(v);
	}
	std::vector<cf_type1>		&m_tc1;
	std::vector<cf_type2>		&m_tc2;
	std::vector<max_fast_int>	&m_ck1;
	std::vector<max_fast_int>	&m_ck2;
	std::vector<term_type1 const *> &m_t1;
	std::vector<term_type2 const *> &m_t2;
	const GenericTruncator		&m_trunc;
	const ArgsTuple			&m_args_tuple;
};


	/// Toolbox for coded series multiplication.
	/**
	 * Intended to be inherited together with piranha::base_series_multiplier. It adds common methods for
	 * series multiplication through Kronecker codification. Requirements:
	 * - input series must have at least one term;
	 * - input series must have at least one argument.
	 */
	template <class Derived, class Series1, class Series2, class OpTuple>
	class coded_multiplier
	{
			// Some static checks.
			p_static_check(Series1::echelon_level == Series2::echelon_level,"");
			p_static_check(boost::tuples::length<OpTuple>::value == Series1::echelon_level + 1,"");
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
			p_static_check((boost::is_same<minmax_type,typename cm_tuple<Series2>::type_minmax>::value),"");
			p_static_check((boost::is_same<value_handler_type,typename cm_tuple<Series2>::type_value_handler>::value),"");
			// Generalised reverse lexicographic comparison.
			class key_revlex_comparison
			{
				public:
					template <class Term>
					bool operator()(const Term *t1, const Term *t2) const
					{
						return key_revlex_comparison_impl<Term>::run(t1,t2);
					}
			};
			enum mult_type
			{
				plain,
				vector,
				hash
			};
			static void trace_mult_type(mult_type type)
			{
				std::string name = "";
				switch (type)
				{
					case plain:
						name = "mult_plain";
						break;
					case vector:
						name = "mult_vector";
						break;
					case hash:
						name = "mult_hash";
				}
				stats::trace_stat(name,std::size_t(0),boost::lambda::_1 + 1);
			}
		public:
			/// Default constructor.
			/**
			 * Initialises viability flag to false and sets up data members for future use. Densities are initialised
			 * to zero and m_mp_gr's and m_fast_gr's sizes are initialised according to the arguments tuple stored in
			 * piranha::base_series_multiplier.
			 */
			coded_multiplier():m_gr_is_viable(false),m_mp_h(mp_integer(0)),m_fast_h(0),
				m_density1(0.),m_density2(0.)
			{
				// NOTE: beware the order of inheritance here, make sure to init base_series_multiplier before,
				//       otherwise m_args_tuple will be uninitialised here.
				// Initialise the member tuples.
				cm_init_vector_tuple<Series1>(m_mp_gr,derived_const_cast->m_args_tuple);
				cm_init_vector_tuple<Series1>(m_fast_gr,derived_const_cast->m_args_tuple);
				cm_init_vector_tuple<Series1>(m_mp_ct,derived_const_cast->m_args_tuple);
				cm_init_vector_tuple<Series1>(m_fast_ct,derived_const_cast->m_args_tuple);
			}
			/// Perform multiplication and place the result into m_retval.
			void perform_multiplication()
			{
				const settings::multiplication_algorithm algo = settings::get_multiplication_algorithm();
				std::vector<typename Series1::term_type> f_terms1;
				std::vector<typename Series2::term_type> f_terms2;
				// If echelon level is more than zero we need to flatten out the series.
				if (Series1::echelon_level) {
					f_terms1 = derived_cast->m_s1.flatten_terms(derived_cast->m_args_tuple);
					f_terms2 = derived_cast->m_s2.flatten_terms(derived_cast->m_args_tuple);
					derived_cast->cache_terms_pointers(f_terms1,f_terms2);
				} else {
					// Cache term pointers.
					derived_cast->cache_terms_pointers(derived_cast->m_s1,derived_cast->m_s2);
				}
				// NOTE: hard coded value of 1000.
				if ((algo == settings::automatic && double(derived_cast->m_terms1.size()) * double(derived_cast->m_terms2.size()) < 1000)
					|| algo == settings::plain)
				{
					derived_cast->perform_plain_multiplication();
					trace_mult_type(plain);
					return;
				}
				// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
				// coefficients and coded keys.
				const typename Derived::truncator_type trunc(derived_cast->m_terms1,derived_cast->m_terms2,derived_cast->m_args_tuple);
				determine_viability();
				if (!m_gr_is_viable) {
					if (algo == settings::vector_coded || algo == settings::hash_coded) {
						piranha_throw(value_error,"coded multiplication requested, but coded representation is infeasible");
					}
					derived_cast->perform_plain_multiplication();
					trace_mult_type(plain);
					return;
				}
				if (trunc.is_effective()) {
					derived_cast->ll_perform_multiplication(trunc);
				} else {
					// Sort input series for better cache usage and multi-threaded implementation.
					std::sort(derived_cast->m_terms1.begin(),derived_cast->m_terms1.end(),key_revlex_comparison());
					std::sort(derived_cast->m_terms2.begin(),derived_cast->m_terms2.end(),key_revlex_comparison());
					derived_cast->ll_perform_multiplication(null_truncator::template get_type<Series1,Series2,typename Derived::args_tuple_type>(
						derived_cast->m_terms1,derived_cast->m_terms2,derived_cast->m_args_tuple
					));
				}
			}
			template <class GenericTruncator>
			void ll_perform_multiplication(const GenericTruncator &trunc)
			{
				const settings::multiplication_algorithm algo = settings::get_multiplication_algorithm();
				typedef typename final_cf<Series1>::type cf_type1;
				typedef typename final_cf<Series2>::type cf_type2;
				// Code terms.
				// NOTE: it is important to code here since at this point we already have sorted input series,
				//       if necessary.
				code_terms();
				// Cache the coefficients.
				std::vector<cf_type1> cf1_cache;
				std::vector<cf_type2> cf2_cache;
				std::insert_iterator<std::vector<cf_type1> > i_it1(cf1_cache,cf1_cache.begin());
				std::insert_iterator<std::vector<cf_type2> > i_it2(cf2_cache,cf2_cache.begin());
				std::transform(derived_cast->m_terms1.begin(),derived_cast->m_terms1.end(),i_it1,final_cf_getter<Series1>());
				std::transform(derived_cast->m_terms2.begin(),derived_cast->m_terms2.end(),i_it2,final_cf_getter<Series2>());
				bool vec_res;
				if ((algo == settings::automatic && is_sparse()) || algo == settings::hash_coded) {
					vec_res = false;
				} else {
					vec_res = derived_cast->perform_vector_coded_multiplication(cf1_cache,cf2_cache,derived_cast->m_terms1,derived_cast->m_terms2,trunc);
				}
				if (!vec_res) {
					if (algo == settings::vector_coded) {
						piranha_throw(value_error,"vector coded multiplication requested, but vector coded representation is infeasible");
					}
					shift_codes();
					derived_cast->perform_hash_coded_multiplication(cf1_cache,cf2_cache,derived_cast->m_terms1,derived_cast->m_terms2,trunc);
					trace_mult_type(hash);
				} else {
					trace_mult_type(vector);
				}
			}
			/// Determine whether the global coded representation is viable or not.
			/**
			 * The m_gr_is_viable flag will be set accordingly after this method is called.
			 */
			void determine_viability()
			{
				// Make sure that the series have at least one term.
				piranha_assert(derived_const_cast->m_terms1.size() > 0 &&
					derived_const_cast->m_terms2.size() > 0);
				// Declare and init the min/max types for the two series.
				minmax_type t1, t2;
				cm_init_vector_tuple<Series1>(t1,derived_const_cast->m_args_tuple);
				cm_init_vector_tuple<Series2>(t2,derived_const_cast->m_args_tuple);
				// Init and test the first series' tuple.
				typedef typename std::vector<typename Series1::term_type const *>::size_type size_type1;
				const size_type1 size1 = derived_const_cast->m_terms1.size();
				// Value handler tuples.
				value_handler_type vh1, vh2;
				cm_minmax<minmax_type>::run_init(*derived_const_cast->m_terms1[0],t1,vh1);
				for (size_type1 i = 1; i < size1; ++i) {
					cm_minmax<minmax_type>::run_test(*derived_const_cast->m_terms1[i],t1,vh1);
				}
				// Init and test the second series' tuple.
				typedef typename std::vector<typename Series2::term_type const *>::size_type size_type2;
				const size_type2 size2 = derived_const_cast->m_terms2.size();
				cm_minmax<minmax_type>::run_init(*derived_const_cast->m_terms2[0],t2,vh2);
				for (size_type2 i = 1; i < size2; ++i) {
					cm_minmax<minmax_type>::run_test(*derived_const_cast->m_terms2[i],t2,vh2);
				}
				// Now compute the global representation in multiprecision.
				cm_global_minmax<OpTuple>::run(t1,vh1,t2,vh2,m_mp_gr);
				// Assign the global value handler tuple.
				piranha_assert(vh1 == vh2);
				m_vh = vh1;
				// Compute the multiprecision coding tuple.
				compute_mp_coding_tuple(m_mp_ct,m_mp_gr);
				// Compute multiprecision codes range.
				tuple_vector_dot(m_mp_gr,m_mp_ct,m_mp_h);
				// To test whether a representation is viable or not, we need to test for the following things:
				// - m_mp_h must be in the max_fast_int range;
				// - m_mp_h's width must be within halft max_fast_int's range (needed for 2*chi shifting).
				// Use lexical cast for max interoperability between numerical types.
				// NOTE: here probably we can reduce greatly the number of memory allocations...
				if (boost::numeric::subset(m_mp_h,boost::numeric::interval<mp_integer>(
					boost::lexical_cast<mp_integer>(boost::integer_traits<max_fast_int>::const_min),
					boost::lexical_cast<mp_integer>(boost::integer_traits<max_fast_int>::const_max))) &&
					boost::numeric::width(m_mp_h) <=
					boost::lexical_cast<mp_integer>(boost::integer_traits<max_fast_int>::const_max) / 2)
				{
					// Mark representation as viable.
					m_gr_is_viable = true;
					// Log viability.
					stats::trace_stat("mult_coded_feasible",std::size_t(0),boost::lambda::_1 + 1);
				} else {
					stats::trace_stat("mult_coded_unfeasible",std::size_t(0),boost::lambda::_1 + 1);
				}
			}
			/// Code terms.
			void code_terms()
			{
				piranha_assert(m_gr_is_viable);
				// Downcast multiprecision to fast representation.
				cm_mp_tuple_downcast(m_mp_gr,m_fast_gr);
				cm_mp_tuple_downcast(m_mp_ct,m_fast_ct);
				m_fast_h.assign(
					boost::lexical_cast<max_fast_int>(m_mp_h.lower()),
					boost::lexical_cast<max_fast_int>(m_mp_h.upper())
				);
				// Build decoding tuple.
				cm_build_decoding_tuple(m_dt,m_fast_gr);
				// Establish if subtraction is requested or not.
				static const bool sub_requested = op_has_sub<OpTuple>::value;
				// Resize codes vectors.
				typedef std::vector<max_fast_int>::size_type size_type;
				const size_type csize1 = boost::numeric_cast<size_type>(derived_const_cast->m_terms1.size()),
					csize2 = boost::numeric_cast<size_type>(derived_const_cast->m_terms2.size());
				m_ckeys1.resize(csize1);
				m_ckeys2a.resize(csize2);
				if (sub_requested) {
					m_ckeys2b.resize(csize2);
				}
				// Now fill in the codes.
				max_fast_int code_a = 0, code_b = 0;
				for (size_type i = 0; i < csize1; ++i) {
					cm_code<OpTuple>(m_fast_ct,*derived_const_cast->m_terms1[i],m_vh,code_a,code_b);
					m_ckeys1[i] = code_a;
				}
				for (size_type i = 0; i < csize2; ++i) {
					cm_code<OpTuple>(m_fast_ct,*derived_const_cast->m_terms2[i],m_vh,code_a,code_b);
					m_ckeys2a[i] = code_a;
					if (sub_requested) {
						m_ckeys2b[i] = code_b;
					}
				}
				// Compute densities.
				const max_fast_int w = boost::numeric::width(m_fast_h) + 1;
				m_density1 = static_cast<double>(csize1) / w;
				m_density2 = static_cast<double>(csize2) / w;
			}
			/// Decode.
			/**
			 * Decode given code into return value term, using final_cf as the coefficient at the end of the echelon recursion.
			 */
			template <class FinalCf>
			void decode(const FinalCf &final_cf, const max_fast_int &code, typename Series1::term_type &term) const
			{
				cm_decode(final_cf,m_dt,m_fast_gr,term,m_vh,code,m_fast_h.lower(),derived_const_cast->m_args_tuple);
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
				piranha_assert(m_gr_is_viable);
				const double max_density = std::max<double>(m_density1,m_density2);
				return (max_density < limit);
			}
			/// Shift codes.
			/**
			 * Move all the the codes so that the minimum code of the representation is 0.
			 */
			void shift_codes()
			{
				piranha_assert(m_gr_is_viable);
				typedef std::vector<max_fast_int>::size_type size_type;
				const size_type size1 = m_ckeys1.size(), size2a = m_ckeys2a.size(), size2b = m_ckeys2b.size();
				const max_fast_int chi = m_fast_h.lower();
				for (size_type i = 0; i < size1; ++i) {
					m_ckeys1[i] -= chi;
					piranha_assert(m_ckeys1[i] >= 0);
				}
				for (size_type i = 0; i < size2a; ++i) {
					m_ckeys2a[i] -= chi;
					piranha_assert(m_ckeys2a[i] >= 0);
				}
				for (size_type i = 0; i < size2b; ++i) {
					m_ckeys2b[i] -= chi;
					piranha_assert(m_ckeys2b[i] >= 0);
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
			boost::numeric::interval<max_fast_int>	m_fast_h;
			/// Codes for the first series.
			std::vector<max_fast_int>		m_ckeys1;
			/// Codes for the second series, plus.
			std::vector<max_fast_int>		m_ckeys2a;
			/// Codes for the second series, minus.
			std::vector<max_fast_int>		m_ckeys2b;
			/// Density of the first series.
			double					m_density1;
			/// Density of the second series.
			double					m_density2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
