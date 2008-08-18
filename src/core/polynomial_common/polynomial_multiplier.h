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
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../p_assert.h"
#include "../settings.h" // For debug.

namespace piranha
{
	/// Series multiplier specifically tuned for Polynomials.
	/**
	 * This multiplier internally will use coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
	class polynomial_multiplier
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
					class term_degree_comparison
					{
						public:
							template <class Term>
							bool operator()(const Term &t1, const Term &t2) const {
								const max_fast_int d1 = t1.m_key.degree(), d2 = t2.m_key.degree();
								if (d1 == d2) {
									return (t1.m_key < t2.m_key);
								} else {
									return (d1 < d2);
								}
							}
					};
				public:
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<get_type> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
							ancestor::base_series_multiplier(s1, s2, retval, args_tuple) {}
					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication() {
						// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
						// coefficients and coded keys.
						const truncator_type trunc(*this);
						coded_ancestor::find_input_min_max();
						calculate_result_min_max();
						coded_ancestor::determine_viability();
						if (trunc.is_effective()) {
							ll_perform_multiplication(trunc);
						} else {
							// We want to sort them this way if we are not truncating, to optimize cache memory usage.
							std::sort(ancestor::m_terms1.begin(),ancestor::m_terms1.end(),term_degree_comparison());
							std::sort(ancestor::m_terms2.begin(),ancestor::m_terms2.end(),term_degree_comparison());
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
							if (coded_ancestor::is_sparse() || !perform_vector_coded_multiplication(trunc)) {
								__PDEBUG(std::cout << "Going for hash coded polynomial multiplication\n");
								perform_hash_coded_multiplication(trunc);
							}
						} else {
							__PDEBUG(std::cout << "Going for plain polynomial multiplication\n");
							ancestor::perform_plain_multiplication(trunc);
						}
					}
					void calculate_result_min_max() {
						std::vector<mpz_class> tmp_vec(6);
						std::pair<typename std::vector<mpz_class>::const_iterator,
							std::vector<mpz_class>::const_iterator> min_max;
						for (size_t i = 0; i < coded_ancestor::m_size; ++i) {
							tmp_vec[0] = coded_ancestor::m_min_max1[i].second;
							tmp_vec[0] += coded_ancestor::m_min_max2[i].second;
							tmp_vec[1] = coded_ancestor::m_min_max1[i].first;
							tmp_vec[1] += coded_ancestor::m_min_max2[i].first;
							tmp_vec[2] = coded_ancestor::m_min_max1[i].first;
							tmp_vec[3] = coded_ancestor::m_min_max2[i].first;
							tmp_vec[4] = coded_ancestor::m_min_max1[i].second;
							tmp_vec[5] = coded_ancestor::m_min_max2[i].second;
							min_max = boost::minmax_element(tmp_vec.begin(), tmp_vec.end());
							coded_ancestor::m_res_min_max[i].first = *(min_max.first);
							coded_ancestor::m_res_min_max[i].second = *(min_max.second);
						}
						__PDEBUG(
							std::cout << "Mult limits are:\n";
						for (size_t i = 0; i < coded_ancestor::m_res_min_max.size(); ++i) {
						std::cout << coded_ancestor::m_res_min_max[i].first << ',' <<
							coded_ancestor::m_res_min_max[i].second << '\n';
						}
						);
					}
					template <class GenericTruncator>
					bool perform_vector_coded_multiplication(const GenericTruncator &trunc) {
						std::vector<cf_type1,std_counting_allocator<cf_type1> > vc;
						// Try to allocate the space for vector coded multiplication.
						// The +1 is needed because we need the number of possible codes between min and max, e.g.:
						// coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
						p_assert(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1 >= 0);
						const size_t n_codes = static_cast<size_t>(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1);
						try {
							vc.resize(n_codes);
						} catch (const std::bad_alloc &) {
							__PDEBUG(std::cout << "Not enough physical memory available for vector coded.\n");
							return false;
						} catch (const out_of_memory &) {
							__PDEBUG(std::cout << "Memory limit reached for vector coded.\n");
							return false;
						}
						__PDEBUG(std::cout << "Going for vector coded polynomial multiplication\n");
						// Define the base pointers for storing the results of multiplication.
						// Please note that even if here it seems like we are going to write outside allocated memory,
						// the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
						cf_type1 *vc_res =  &vc[0] - coded_ancestor::m_h_min;
						// Perform multiplication.
						for (size_t i = 0; i < ancestor::m_size1; ++i) {
							const max_fast_int index1 = coded_ancestor::m_ckeys1[i];
							for (size_t j = 0; j < ancestor::m_size2; ++j) {
								// Calculate index of the result.
								const max_fast_int res_index = index1 + coded_ancestor::m_ckeys2[j];
								if (trunc.skip(ancestor::m_terms1[i], ancestor::m_terms2[j])) {
									break;
								}
								if (trunc.accept(res_index)) {
									vc_res[res_index].addmul(ancestor::m_terms1[i].m_cf, ancestor::m_terms2[j].m_cf,
															 ancestor::m_args_tuple);
								}
							}
						}
						__PDEBUG(std::cout << "Done multiplying\n");
						// Decode and insert the results into return value.
						term_type1 tmp_term;
						for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i) {
							// Take a shortcut and check for ignorability of the coefficient here.
							// This way we avoid decodification, and all the series term insertion yadda-yadda.
							if (!vc_res[i].is_ignorable(ancestor::m_args_tuple)) {
								tmp_term.m_cf = vc_res[i];
								coded_ancestor::decode(tmp_term.m_key, i);
								if (!tmp_term.is_canonical(ancestor::m_args_tuple)) {
									tmp_term.canonicalise(ancestor::m_args_tuple);
								}
								ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done polynomial vector coded\n");
						return true;
					}
					template <class GenericTruncator>
					void perform_hash_coded_multiplication(const GenericTruncator &trunc) {
						typedef coded_series_cuckoo_hash_table<cf_type1, max_fast_int, std_counting_allocator<char> > csht;
						typedef typename csht::term_type cterm;
						typedef typename csht::iterator c_iterator;
						// Let's find a sensible size hint.
						const size_t size_hint = std::max<double>(this->m_density1,this->m_density2) * this->m_h_tot;
						csht cms(size_hint);
						cterm tmp_cterm;
						for (size_t i = 0; i < ancestor::m_size1; ++i) {
							const max_fast_int key1 = coded_ancestor::m_ckeys1[i];
							for (size_t j = 0; j < ancestor::m_size2; ++j) {
								if (trunc.skip(ancestor::m_terms1[i], ancestor::m_terms2[j])) {
									break;
								}
								const max_fast_int new_key = key1 + coded_ancestor::m_ckeys2[j];
								switch (trunc.accept(new_key)) {
								case true: {
									c_iterator it = cms.find(new_key);
									switch (it == cms.end()) {
									case true: {
										// Create new temporary term from old cf and new key.
										tmp_cterm.m_cf = ancestor::m_terms1[i].m_cf;
										tmp_cterm.m_ckey = new_key;
										// Multiply the old term by the second term.
										tmp_cterm.m_cf.mult_by(ancestor::m_terms2[j].m_cf, ancestor::m_args_tuple);
										cms.insert(tmp_cterm);
										break;
									}
									case false:
										it->m_cf.addmul(ancestor::m_terms1[i].m_cf, ancestor::m_terms2[j].m_cf,
														ancestor::m_args_tuple);
									}
								}
								case false:
									;
								}
							}
						}
						__PDEBUG(std::cout << "Done polynomial hash coded multiplying\n");
						// Decode and insert into retval.
						ancestor::m_retval.rehash(cms.size());
						term_type1 tmp_term;
						const c_iterator c_it_f = cms.end();
						for (c_iterator c_it = cms.begin(); c_it != c_it_f; ++c_it) {
							tmp_term.m_cf = c_it->m_cf;
							coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
							if (!tmp_term.is_canonical(ancestor::m_args_tuple)) {
								tmp_term.canonicalise(ancestor::m_args_tuple);
							}
							ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
						}
						__PDEBUG(std::cout << "Done polynomial hash coded\n");
					}
				public:
					// Temporary key used for the decodification in the truncator.
					// It is mutable because it is used as temporary decodification area.
					mutable typename term_type1::key_type m_tmp_key;
			};
	};
}

#endif
