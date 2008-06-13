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

#ifndef PIRANHA_POISSON_SERIES_MULTIPLIER_H
#define PIRANHA_POISSON_SERIES_MULTIPLIER_H

#include <boost/algorithm/minmax_element.hpp> // To calculate limits of multiplication.
#include <boost/static_assert.hpp>
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_series_multiplier.h"
#include "../base_classes/coded_series_hash_table.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../settings.h" // For debug.

namespace piranha
{
	/// Series multiplier specifically tuned for Poisson series.
	/**
	 * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
	template <class Series1, class Series2, class ArgsTuple, template <class> class Truncator>
	class poisson_series_multiplier:
				public base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
				poisson_series_multiplier<Series1, Series2, ArgsTuple, Truncator> > ,
				public coded_series_multiplier<poisson_series_multiplier<Series1, Series2, ArgsTuple, Truncator> >
	{
			typedef base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
			poisson_series_multiplier<Series1, Series2, ArgsTuple, Truncator> > ancestor;
			typedef coded_series_multiplier<poisson_series_multiplier<Series1, Series2, ArgsTuple, Truncator> > coded_ancestor;
			friend class coded_series_multiplier<poisson_series_multiplier<Series1, Series2, ArgsTuple, Truncator> >;
			friend class Truncator<poisson_series_multiplier<Series1, Series2, ArgsTuple, Truncator> >;
			typedef typename Series1::const_sorted_iterator const_iterator1;
			typedef typename Series2::const_sorted_iterator const_iterator2;
			typedef typename Series1::sorted_iterator iterator1;
			typedef typename Series2::sorted_iterator iterator2;
			typedef Series1 series_type1;
			typedef Series2 series_type2;
			typedef typename ancestor::term_type1 term_type1;
			typedef typename ancestor::term_type2 term_type2;
			typedef typename term_type1::cf_type cf_type1;
			typedef typename term_type2::cf_type cf_type2;
			typedef typename term_type1::key_type key_type;
		public:
			typedef typename ancestor::truncator_type truncator_type;
			poisson_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
					ancestor::base_series_multiplier(s1, s2, retval, args_tuple),
					m_cfs1(ancestor::m_size1), m_cfs2(ancestor::m_size2),
					m_flavours1(ancestor::m_size1), m_flavours2(ancestor::m_size2) {}
			/// Perform multiplication and place the result into m_retval.
			void perform_multiplication() {
				coded_ancestor::find_input_min_max();
				calculate_result_min_max();
				coded_ancestor::determine_viability();
				if (coded_ancestor::m_cr_is_viable) {
					// Here we should be ok, since we know that the two sizes are greater than zero and even
					// if we divide by zero we should get Inf, which is fine for our purposes.
					const double density = ((double)ancestor::m_size1 * ancestor::m_size2 * 2) /
										   (coded_ancestor::m_h_max - coded_ancestor::m_h_min);
					__PDEBUG(std::cout << "Density: " << density << '\n');
					coded_ancestor::code_keys();
					cache_flavours();
					cache_coefficients();
					if (density < 1E-1 or !perform_vector_coded_multiplication()) {
						__PDEBUG(if (density < 1E-1) std::cout << "Low density\n");
						__PDEBUG(std::cout << "Going for hash coded Poisson series multiplication\n");
						perform_hash_coded_multiplication();
					}
				} else {
					__PDEBUG(std::cout << "Going for plain Poisson series multiplication\n");
					ancestor::perform_plain_multiplication();
				}
			}
		private:
			void calculate_result_min_max() {
				// TODO: optimize here the usage of mpz classes.
				std::vector<mpz_class> tmp_vec(8);
				std::pair<typename std::vector<mpz_class>::const_iterator, std::vector<mpz_class>::const_iterator> min_max;
				for (size_t i = 0; i < coded_ancestor::m_size; ++i) {
					tmp_vec[0] = coded_ancestor::m_min_max1[i].second;
					tmp_vec[0] += coded_ancestor::m_min_max2[i].second;
					tmp_vec[1] = coded_ancestor::m_min_max1[i].first;
					tmp_vec[1] += coded_ancestor::m_min_max2[i].first;
					tmp_vec[2] = coded_ancestor::m_min_max1[i].second;
					tmp_vec[2] -= coded_ancestor::m_min_max2[i].first;
					tmp_vec[3] = coded_ancestor::m_min_max1[i].first;
					tmp_vec[3] -= coded_ancestor::m_min_max2[i].second;
					tmp_vec[4] = coded_ancestor::m_min_max1[i].first;
					tmp_vec[5] = coded_ancestor::m_min_max2[i].first;
					tmp_vec[6] = coded_ancestor::m_min_max1[i].second;
					tmp_vec[7] = coded_ancestor::m_min_max2[i].second;
					min_max = boost::minmax_element(tmp_vec.begin(), tmp_vec.end());
					coded_ancestor::m_res_min_max[i].first = *(min_max.first);
					coded_ancestor::m_res_min_max[i].second = *(min_max.second);
				}
				__PDEBUG(
					std::cout << "Mult limits are:\n";
				for (size_t i = 0; i < coded_ancestor::m_res_min_max.size(); ++i) {
				std::cout << coded_ancestor::m_res_min_max[i].first << ',' << coded_ancestor::m_res_min_max[i].second << '\n';
				}
				);
			}
			/// Store flavours of the series into own vectors.
			void cache_flavours() {
				size_t i;
				for (i = 0; i < ancestor::m_size1; ++i) {
					m_flavours1[i] = ancestor::m_terms1[i]->m_key.flavour();
				}
				for (i = 0; i < ancestor::m_size2; ++i) {
					m_flavours2[i] = ancestor::m_terms2[i]->m_key.flavour();
				}
			}
			void cache_coefficients() {
				size_t i;
				for (i = 0; i < ancestor::m_size1; ++i) {
					m_cfs1[i] = ancestor::m_terms1[i]->m_cf;
				}
				for (i = 0; i < ancestor::m_size2; ++i) {
					m_cfs2[i] = ancestor::m_terms2[i]->m_cf;
				}
			}
			bool perform_vector_coded_multiplication() {
				cf_type1 *p_vc_res_cos(0), *p_vc_res_sin(0);
				// Try to allocate the space for vector coded multiplication. We need two arrays of results,
				// one for cosines, one for sines.
				// The +1 is needed because we need the number of possible codes between min and max, e.g.:
				// coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
				const size_t n_codes = (size_t)(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1);
				try {
					p_vc_res_cos = (cf_type1 *)piranha_malloc(sizeof(cf_type1) * n_codes * 2);
					p_vc_res_sin = p_vc_res_cos + n_codes;
					// Reset memory area. Use positional new so that if cf is a class with non-trivial ctors,
					// we are sure it will be initialized properly. We want to make sure the coefficients are initialized
					// to zero in order to accumulate Poisson terms during multiplication.
					for (size_t i = 0; i < n_codes; ++i) {
						::new(p_vc_res_cos + i) cf_type1((max_fast_int)0, ancestor::m_args_tuple);
						::new(p_vc_res_sin + i) cf_type1((max_fast_int)0, ancestor::m_args_tuple);
					}
				} catch (const std::bad_alloc &) {
					piranha_free(p_vc_res_cos);
					return false;
				}
				__PDEBUG(std::cout << "Going for vector coded Poisson series multiplication\n");
				// Define the base pointers for storing the results of multiplication.
				// Please note that even if here it seems like we are going to write outside allocated memory,
				// the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
				cf_type1 *vc_res_cos =  p_vc_res_cos - coded_ancestor::m_h_min, *vc_res_sin = p_vc_res_sin - coded_ancestor::m_h_min;
				cf_type1 tmp_cf;
				// Perform multiplication.
				for (size_t i = 0; i < ancestor::m_size1; ++i) {
					for (size_t j = 0; j < ancestor::m_size2; ++j) {
						if (ancestor::m_trunc.skip(ancestor::m_terms1[i], ancestor::m_terms2[j])) {
							break;
						}
						// TODO: Does it make sense here to define a method for coefficients like:
						// mult_by_and_insert_into<bool Sign>(cf2,retval,m_args_tuple)
						// so that we can avoid copying stuff around here and elsewhere?
						tmp_cf = m_cfs1[i];
						tmp_cf.mult_by(m_cfs2[j], ancestor::m_args_tuple);
						tmp_cf.divide_by((max_fast_int)2, ancestor::m_args_tuple);
						const max_fast_int index_plus = coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j],
														index_minus = coded_ancestor::m_ckeys1[i] - coded_ancestor::m_ckeys2[j];
						switch (m_flavours1[i] == m_flavours2[j]) {
						case true:
							switch (m_flavours1[i]) {
							case true:
								vc_res_cos[index_minus].add(tmp_cf, ancestor::m_args_tuple);
								vc_res_cos[index_plus].add(tmp_cf, ancestor::m_args_tuple);
								break;
							case false:
								vc_res_cos[index_minus].add(tmp_cf, ancestor::m_args_tuple);
								vc_res_cos[index_plus].subtract(tmp_cf, ancestor::m_args_tuple);
							}
							break;
						case false:
							switch (m_flavours1[i]) {
							case true:
								vc_res_sin[index_minus].subtract(tmp_cf, ancestor::m_args_tuple);
								vc_res_sin[index_plus].add(tmp_cf, ancestor::m_args_tuple);
								break;
							case false:
								vc_res_sin[index_minus].add(tmp_cf, ancestor::m_args_tuple);
								vc_res_sin[index_plus].add(tmp_cf, ancestor::m_args_tuple);
							}
						}
					}
				}
				__PDEBUG(std::cout << "Done multiplying\n");
				// Decode and insert the results into return value.
				term_type1 tmp_term;
				iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
				for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i) {
					// Take a shortcut and check for ignorability of the coefficient here.
					// This way we avoid decodification, and all the series term insertion yadda-yadda.
					switch (likely(vc_res_cos[i].is_ignorable(ancestor::m_args_tuple))) {
					case true:
						break;
					case false:
						tmp_term.m_cf = vc_res_cos[i];
						coded_ancestor::decode(tmp_term.m_key, i);
						tmp_term.m_key.flavour() = true;
						it_hint = ancestor::m_retval.insert(tmp_term, it_hint, ancestor::m_args_tuple);
					}
				}
				for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i) {
					switch (likely(vc_res_sin[i].is_ignorable(ancestor::m_args_tuple))) {
					case true:
						break;
					case false:
						tmp_term.m_cf = vc_res_sin[i];
						coded_ancestor::decode(tmp_term.m_key, i);
						tmp_term.m_key.flavour() = false;
						it_hint = ancestor::m_retval.insert(tmp_term, it_hint, ancestor::m_args_tuple);
					}
				}
				// Call dtors for the coefficients in the allocated space.
				// This is necessary for non-trivial coefficients.
				for (size_t i = 0; i < n_codes; ++i) {
					(p_vc_res_cos + i)->~cf_type1();
					(p_vc_res_sin + i)->~cf_type1();
				}
				// Free the allocated space.
				piranha_free(p_vc_res_cos);
				__PDEBUG(std::cout << "Done Poisson series vector coded\n");
				return true;
			}
			void perform_hash_coded_multiplication() {
				typedef coded_series_hash_table<cf_type1, max_fast_int> csht;
				typedef typename csht::term_type cterm;
				typedef typename csht::iterator c_iterator;
				csht cms_cos, cms_sin;
				for (size_t i = 0; i < ancestor::m_size1; ++i) {
					for (size_t j = 0; j < ancestor::m_size2; ++j) {
						if (ancestor::m_trunc.skip(ancestor::m_terms1[i], ancestor::m_terms2[j])) {
							break;
						}
						// TODO: here (and elsewhere, likely), we can avoid an extra copy by working with keys and cfs instead of terms,
						// generating only one coefficient and change its sign later if needed - after insertion.
						cterm tmp_term1(m_cfs1[i], coded_ancestor::m_ckeys1[i]);
						// Handle the coefficient, with positive signs for now.
						tmp_term1.m_cf.mult_by(m_cfs2[j], ancestor::m_args_tuple);
						tmp_term1.m_cf.divide_by((max_fast_int)2, ancestor::m_args_tuple);
						tmp_term1.m_ckey -= coded_ancestor::m_ckeys2[j];
						// Create the second term, using the first one's coefficient and the appropriate code.
						cterm tmp_term2(tmp_term1.m_cf, coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j]);
						// Now fix flavours and coefficient signs.
						switch (m_flavours1[i] == m_flavours2[j]) {
						case true: {
							switch (m_flavours1[i]) {
							case true:
								break;
							case false:
								tmp_term2.m_cf.invert_sign(ancestor::m_args_tuple);
							}
							// Insert into cosine container.
							c_iterator it = cms_cos.find(tmp_term1.m_ckey);
							switch (it == cms_cos.end()) {
							case true:
								cms_cos.insert(tmp_term1);
								break;
							case false:
								it->m_cf.add(tmp_term1.m_cf, ancestor::m_args_tuple);
							}
							it = cms_cos.find(tmp_term2.m_ckey);
							switch (it == cms_cos.end()) {
							case true:
								cms_cos.insert(tmp_term2);
								break;
							case false:
								it->m_cf.add(tmp_term2.m_cf, ancestor::m_args_tuple);
							}
							break;
						}
						case false: {
							switch (m_flavours1[i]) {
							case true:
								tmp_term1.m_cf.invert_sign(ancestor::m_args_tuple);
								break;
							case false:
								;
							}
							// Insert into sine container.
							c_iterator it = cms_sin.find(tmp_term1.m_ckey);
							switch (it == cms_sin.end()) {
							case true:
								cms_sin.insert(tmp_term1);
								break;
							case false:
								it->m_cf.add(tmp_term1.m_cf, ancestor::m_args_tuple);
							}
							it = cms_sin.find(tmp_term2.m_ckey);
							switch (it == cms_sin.end()) {
							case true:
								cms_sin.insert(tmp_term2);
								break;
							case false:
								it->m_cf.add(tmp_term2.m_cf, ancestor::m_args_tuple);
							}
						}
						}
					}
				}
				__PDEBUG(std::cout << "Done Poisson series hash coded multiplying\n");
				term_type1 tmp_term;
				iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
				{
					const c_iterator c_it_f = cms_cos.end();
					for (c_iterator c_it = cms_cos.begin(); c_it != c_it_f; ++c_it) {
						tmp_term.m_cf = c_it->m_cf;
						coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
						tmp_term.m_key.flavour() = true;
						it_hint = ancestor::m_retval.insert(tmp_term, it_hint, ancestor::m_args_tuple);
					}
				}
				{
					const c_iterator c_it_f = cms_sin.end();
					for (c_iterator c_it = cms_sin.begin(); c_it != c_it_f; ++c_it) {
						tmp_term.m_cf = c_it->m_cf;
						coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
						tmp_term.m_key.flavour() = false;
						it_hint = ancestor::m_retval.insert(tmp_term, it_hint, ancestor::m_args_tuple);
					}
				}
				__PDEBUG(std::cout << "Done Poisson series hash coded\n");
			}
		private:
			std::vector<cf_type1>	m_cfs1;
			std::vector<cf_type2>	m_cfs2;
			// For Poisson series we also need flavours.
			std::vector<char>		m_flavours1;
			std::vector<char>		m_flavours2;
			// Just making sure, eh...
			BOOST_STATIC_ASSERT(sizeof(char) == sizeof(bool));
	};
}

#endif
