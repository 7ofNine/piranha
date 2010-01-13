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
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp> // We assert equality between vh tuples below.
#include <boost/type_traits/is_same.hpp>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"
#include "coded_multiplier_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
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
		public:
			template <class Cf, class Ckey>
			struct coded_term_type {
				coded_term_type():m_cf(),m_ckey() {}
				coded_term_type(const Cf &cf, const Ckey &ckey):m_cf(cf),m_ckey(ckey) {}
				bool operator==(const coded_term_type &t) const
				{
					return (m_ckey == t.m_ckey);
				}
				std::size_t hash_value() const
				{
					return boost::hash<Ckey>()(m_ckey);
				}
				mutable Cf	m_cf;
				Ckey		m_ckey;
			};
		protected:
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
				// - m_mp_h's width must be in the max_fast_int range.
				// Use lexical cast for max interoperability between numerical types.
				// NOTE: here probably we can reduce greatly the number of memory allocations...
				if (boost::numeric::subset(m_mp_h,boost::numeric::interval<mp_integer>(
					boost::lexical_cast<mp_integer>(boost::integer_traits<max_fast_int>::const_min),
					boost::lexical_cast<mp_integer>(boost::integer_traits<max_fast_int>::const_max))) &&
					boost::numeric::width(m_mp_h) <=
					boost::lexical_cast<mp_integer>(boost::integer_traits<max_fast_int>::const_max)) {
					// Mark representation as viable.
					m_gr_is_viable = true;
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
