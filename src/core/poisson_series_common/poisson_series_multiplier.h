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
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_series_multiplier.h"
#include "../base_classes/coded_series_cuckoo_hash_table.h"
#include "../base_classes/null_truncator.h"
#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../p_assert.h"
#include "../settings.h" // For debug.

namespace piranha
{
	/// Series multiplier specifically tuned for Poisson series.
	/**
	 * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
	// TODO: isn't truncator's "accept" used here??
	class poisson_series_multiplier
	{
		public:
			template <class Series1, class Series2, class ArgsTuple, class Truncator>
			class get_type:
						public base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
						get_type<Series1, Series2, ArgsTuple, Truncator> > ,
						public coded_series_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator> >
			{
					typedef base_series_multiplier < Series1, Series2, ArgsTuple, Truncator,
					get_type<Series1, Series2, ArgsTuple, Truncator> > ancestor;
					typedef coded_series_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator> > coded_ancestor;
					friend class coded_series_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator> >;
					typedef typename Series1::const_iterator const_iterator1;
					typedef typename Series2::const_iterator const_iterator2;
					typedef typename ancestor::term_type1 term_type1;
					typedef typename ancestor::term_type2 term_type2;
					typedef typename term_type1::cf_type cf_type1;
					typedef typename term_type2::cf_type cf_type2;
					typedef typename term_type1::key_type key_type;
				public:
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<get_type> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
							ancestor::base_series_multiplier(s1, s2, retval, args_tuple),
							m_flavours1(ancestor::m_size1), m_flavours2(ancestor::m_size2) {}
					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication() {
						// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
						// coefficients and coded keys.
						const truncator_type trunc(*this);
						coded_ancestor::find_input_min_max();
						calculate_result_min_max();
						coded_ancestor::determine_viability();
						// Use the selected truncator only if it really truncates, otherwise use the
						// null truncator.
						if (trunc.is_effective()) {
							ll_perform_multiplication(trunc);
						} else {
							ll_perform_multiplication(null_truncator::template get_type<get_type>(*this));
						}
					}
				private:
					template <class GenericTruncator>
					void ll_perform_multiplication(const GenericTruncator &trunc) {
						if (ancestor::m_terms1.size() < 10 && ancestor::m_terms2.size() < 10) {
							__PDEBUG(std::cout << "Small series, going for plain polynomial multiplication\n");
							ancestor::perform_plain_multiplication(trunc);
						} else if (coded_ancestor::m_cr_is_viable) {
							coded_ancestor::code_keys();
							cache_flavours();
							if (coded_ancestor::is_sparse() || !perform_vector_coded_multiplication(trunc)) {
								__PDEBUG(std::cout << "Going for hash coded Poisson series multiplication\n");
								perform_hash_coded_multiplication(trunc);
							}
						} else {
							__PDEBUG(std::cout << "Going for plain Poisson series multiplication\n");
							ancestor::perform_plain_multiplication(trunc);
						}
					}
					void calculate_result_min_max() {
						std::vector<mpz_class> tmp_vec(8);
						std::pair<typename std::vector<mpz_class>::const_iterator,
							std::vector<mpz_class>::const_iterator> min_max;
						const size_t size = coded_ancestor::m_size;
						for (size_t i = 0; i < size; ++i) {
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
							m_flavours1[i] = ancestor::m_terms1[i].m_key.flavour();
						}
						for (i = 0; i < ancestor::m_size2; ++i) {
							m_flavours2[i] = ancestor::m_terms2[i].m_key.flavour();
						}
					}
					template <class GenericTruncator>
					bool perform_vector_coded_multiplication(const GenericTruncator &trunc) {
						std::vector<cf_type1,std_counting_allocator<cf_type1> > vc;
						// Try to allocate the space for vector coded multiplication. We need two arrays of results,
						// one for cosines, one for sines.
						// The +1 is needed because we need the number of possible codes between min and max, e.g.:
						// coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
						p_assert(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1 >= 0);
						const size_t n_codes = static_cast<size_t>(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1);
						try {
							vc.resize(n_codes << 1);
						} catch (const std::bad_alloc &) {
							__PDEBUG(std::cout << "Not enough physical memory available for vector coded.\n");
							return false;
						} catch (const out_of_memory &) {
							__PDEBUG(std::cout << "Memory limit reached for vector coded.\n");
							return false;
						}
						__PDEBUG(std::cout << "Going for vector coded Poisson series multiplication\n");
						// Define the base pointers for storing the results of multiplication.
						// Please note that even if here it seems like we are going to write outside allocated memory,
						// the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
						cf_type1 *vc_res_cos =  &vc[0] - coded_ancestor::m_h_min,
							*vc_res_sin = &vc[0] + n_codes - coded_ancestor::m_h_min;
						cf_type1 tmp_cf;
						// Perform multiplication.
						// TODO: for better cache behaviour and to reduce branching, maybe we can split up input series
						// into cosine / sine parts and multiply them separately. Also, we can do separately
						// index minus and index plus.
						for (size_t i = 0; i < ancestor::m_size1; ++i) {
							for (size_t j = 0; j < ancestor::m_size2; ++j) {
								if (trunc.skip(ancestor::m_terms1[i], ancestor::m_terms2[j])) {
									break;
								}
								// TODO: Does it make sense here to define a method for coefficients like:
								// mult_by_and_insert_into<bool Sign>(cf2,retval,m_args_tuple)
								// so that we can avoid copying stuff around here and elsewhere?
								tmp_cf = ancestor::m_terms1[i].m_cf;
								tmp_cf.mult_by(ancestor::m_terms2[j].m_cf, ancestor::m_args_tuple);
								tmp_cf.divide_by(static_cast<max_fast_int>(2), ancestor::m_args_tuple);
								const max_fast_int index_plus = coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j],
																index_minus = coded_ancestor::m_ckeys1[i] - coded_ancestor::m_ckeys2[j];
								if (m_flavours1[i] == m_flavours2[j]) {
									if (m_flavours1[i]) {
										vc_res_cos[index_minus].add(tmp_cf, ancestor::m_args_tuple);
										vc_res_cos[index_plus].add(tmp_cf, ancestor::m_args_tuple);
									} else {
										vc_res_cos[index_minus].add(tmp_cf, ancestor::m_args_tuple);
										vc_res_cos[index_plus].subtract(tmp_cf, ancestor::m_args_tuple);
									}
								} else {
									if (m_flavours1[i]) {
										vc_res_sin[index_minus].subtract(tmp_cf, ancestor::m_args_tuple);
										vc_res_sin[index_plus].add(tmp_cf, ancestor::m_args_tuple);
									} else {
										vc_res_sin[index_minus].add(tmp_cf, ancestor::m_args_tuple);
										vc_res_sin[index_plus].add(tmp_cf, ancestor::m_args_tuple);
									}
								}
							}
						}
						__PDEBUG(std::cout << "Done multiplying\n");
						// Decode and insert the results into return value.
						term_type1 tmp_term;
						for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i) {
							// Take a shortcut and check for ignorability of the coefficient here.
							// This way we avoid decodification, and all the series term insertion yadda-yadda.
							if (!vc_res_cos[i].is_ignorable(ancestor::m_args_tuple)) {
								tmp_term.m_cf = vc_res_cos[i];
								coded_ancestor::decode(tmp_term.m_key, i);
								tmp_term.m_key.flavour() = true;
								// Canonicalise in-place, so that we don't need to make further copies in the
								// main insertion function.
								if (!tmp_term.is_canonical(ancestor::m_args_tuple)) {
									tmp_term.canonicalise(ancestor::m_args_tuple);
								}
								ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
							}
						}
						for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i) {
							if (!vc_res_sin[i].is_ignorable(ancestor::m_args_tuple)) {
								tmp_term.m_cf = vc_res_sin[i];
								coded_ancestor::decode(tmp_term.m_key, i);
								tmp_term.m_key.flavour() = false;
								if (!tmp_term.is_canonical(ancestor::m_args_tuple)) {
									tmp_term.canonicalise(ancestor::m_args_tuple);
								}
								ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done Poisson series vector coded\n");
						return true;
					}
					template <class GenericTruncator>
					void perform_hash_coded_multiplication(const GenericTruncator &trunc) {
						typedef coded_series_cuckoo_hash_table<cf_type1, max_fast_int, std_counting_allocator<char> > csht;
						typedef typename csht::term_type cterm;
						typedef typename csht::iterator c_iterator;
						csht cms_cos, cms_sin;
						for (size_t i = 0; i < ancestor::m_size1; ++i) {
							for (size_t j = 0; j < ancestor::m_size2; ++j) {
								if (trunc.skip(ancestor::m_terms1[i], ancestor::m_terms2[j])) {
									break;
								}
								// TODO: here (and elsewhere, likely), we can avoid an extra copy by working with keys and cfs instead of terms,
								// generating only one coefficient and change its sign later if needed - after insertion.
								cterm tmp_term1(ancestor::m_terms1[i].m_cf, coded_ancestor::m_ckeys1[i]);
								// Handle the coefficient, with positive signs for now.
								tmp_term1.m_cf.mult_by(ancestor::m_terms2[j].m_cf, ancestor::m_args_tuple);
								tmp_term1.m_cf.divide_by(static_cast<max_fast_int>(2), ancestor::m_args_tuple);
								tmp_term1.m_ckey -= coded_ancestor::m_ckeys2[j];
								// Create the second term, using the first one's coefficient and the appropriate code.
								cterm tmp_term2(tmp_term1.m_cf, coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j]);
								// Now fix flavours and coefficient signs.
								if (m_flavours1[i] == m_flavours2[j]) {
									if (!m_flavours1[i]) {
										tmp_term2.m_cf.invert_sign(ancestor::m_args_tuple);
									}
									// Insert into cosine container.
									c_iterator it = cms_cos.find(tmp_term1.m_ckey);
									if (it == cms_cos.end()) {
										cms_cos.insert(tmp_term1);
									} else {
										it->m_cf.add(tmp_term1.m_cf, ancestor::m_args_tuple);
									}
									it = cms_cos.find(tmp_term2.m_ckey);
									if (it == cms_cos.end()) {
										cms_cos.insert(tmp_term2);
									} else {
										it->m_cf.add(tmp_term2.m_cf, ancestor::m_args_tuple);
									}
								} else {
									if (m_flavours1[i]) {
										tmp_term1.m_cf.invert_sign(ancestor::m_args_tuple);
									}
									// Insert into sine container.
									c_iterator it = cms_sin.find(tmp_term1.m_ckey);
									if (it == cms_sin.end()) {
										cms_sin.insert(tmp_term1);
									} else {
										it->m_cf.add(tmp_term1.m_cf, ancestor::m_args_tuple);
									}
									it = cms_sin.find(tmp_term2.m_ckey);
									if (it == cms_sin.end()) {
										cms_sin.insert(tmp_term2);
									} else {
										it->m_cf.add(tmp_term2.m_cf, ancestor::m_args_tuple);
									}
								}
							}
						}
						__PDEBUG(std::cout << "Done Poisson series hash coded multiplying\n");
						ancestor::m_retval.rehash(cms_cos.size() + cms_sin.size());
						term_type1 tmp_term;
						{
							const c_iterator c_it_f = cms_cos.end();
							for (c_iterator c_it = cms_cos.begin(); c_it != c_it_f; ++c_it) {
								tmp_term.m_cf = c_it->m_cf;
								coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
								tmp_term.m_key.flavour() = true;
								if (!tmp_term.is_canonical(ancestor::m_args_tuple)) {
									tmp_term.canonicalise(ancestor::m_args_tuple);
								}
								ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
							}
						}
						{
							const c_iterator c_it_f = cms_sin.end();
							for (c_iterator c_it = cms_sin.begin(); c_it != c_it_f; ++c_it) {
								tmp_term.m_cf = c_it->m_cf;
								coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
								tmp_term.m_key.flavour() = false;
								if (!tmp_term.is_canonical(ancestor::m_args_tuple)) {
									tmp_term.canonicalise(ancestor::m_args_tuple);
								}
								ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done Poisson series hash coded\n");
					}
				private:
					// For Poisson series we also need flavours.
					std::vector<char>		m_flavours1;
					std::vector<char>		m_flavours2;
			};
	};
}

#endif
