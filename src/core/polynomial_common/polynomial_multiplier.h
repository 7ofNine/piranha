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

#include <boost/algorithm/minmax_element.hpp> // To calculate limits of multiplication.
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
					friend class Truncator::template get_type<get_type>::type;
					typedef typename Series1::template const_iterator<0>::type const_iterator1;
					typedef typename Series2::template const_iterator<0>::type const_iterator2;
					typedef Series1 series_type1;
					typedef Series2 series_type2;
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
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<get_type> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
							ancestor::base_series_multiplier(s1, s2, retval, args_tuple) {
// 						std::sort(ancestor::m_terms1.begin(),ancestor::m_terms1.end(),term_degree_comparison());
// 						std::sort(ancestor::m_terms2.begin(),ancestor::m_terms2.end(),term_degree_comparison());
					}
					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication() {
						coded_ancestor::find_input_min_max();
						calculate_result_min_max();
						coded_ancestor::determine_viability();
						// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
						// coefficients and coded keys.
						const truncator_type trunc(*this);
						if (coded_ancestor::m_cr_is_viable) {
							// Here we should be ok, since we know that the two sizes are greater than zero and even
							// if we divide by zero we should get Inf, which is fine for our purposes.
							const double density = ((double)ancestor::m_size1 * ancestor::m_size2) /
												   (coded_ancestor::m_h_max - coded_ancestor::m_h_min);
							__PDEBUG(std::cout << "Density: " << density << '\n');
							coded_ancestor::code_keys();
							if (density < 1E-1 || !perform_vector_coded_multiplication(trunc)) {
								__PDEBUG(if (density < 1E-1) std::cout << "Low density\n");
								__PDEBUG(std::cout << "Going for hash coded polynomial multiplication\n");
								perform_hash_coded_multiplication(trunc);
							}
						} else {
							__PDEBUG(std::cout << "Going for plain polynomial multiplication\n");
							ancestor::perform_plain_multiplication(trunc);
						}
					}
				private:
					void calculate_result_min_max() {
						std::vector<mpz_class> tmp_vec(6);
						std::pair<typename std::vector<mpz_class>::const_iterator, std::vector<mpz_class>::const_iterator> min_max;
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
						std::cout << coded_ancestor::m_res_min_max[i].first << ',' << coded_ancestor::m_res_min_max[i].second << '\n';
						}
						);
					}
					template <class GenericTruncator>
					bool perform_vector_coded_multiplication(const GenericTruncator &trunc) {
						cf_type1 *p_vc_res(0);
						// Try to allocate the space for vector coded multiplication.
						// The +1 is needed because we need the number of possible codes between min and max, e.g.:
						// coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
						const size_t n_codes = (size_t)(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1);
						try {
							p_vc_res = (cf_type1 *)piranha_malloc(sizeof(cf_type1) * n_codes);
							// Reset memory area. Use positional new so that if cf is a class with non-trivial ctors,
							// we are sure it will be initialized properly. We want to make sure the coefficients are initialized
							// to zero in order to accumulate monomials during multiplication.
							for (size_t i = 0; i < n_codes; ++i) {
								::new(p_vc_res + i) cf_type1((max_fast_int)0, ancestor::m_args_tuple);
							}
						} catch (const std::bad_alloc &) {
							piranha_free(p_vc_res);
							return false;
						}
						__PDEBUG(std::cout << "Going for vector coded polynomial multiplication\n");
						// Define the base pointers for storing the results of multiplication.
						// Please note that even if here it seems like we are going to write outside allocated memory,
						// the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
						cf_type1 *vc_res =  p_vc_res - coded_ancestor::m_h_min;
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
							switch (likely(vc_res[i].is_ignorable(ancestor::m_args_tuple))) {
							case true:
								break;
							case false:
								tmp_term.m_cf = vc_res[i];
								coded_ancestor::decode(tmp_term.m_key, i);
								ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
							}
						}
						// Call dtors for the coefficients in the allocated space.
						// This is necessary for non-trivial coefficients.
						for (size_t i = 0; i < n_codes; ++i) {
							(p_vc_res + i)->~cf_type1();
						}
						// Free the allocated space.
						piranha_free(p_vc_res);
						__PDEBUG(std::cout << "Done polynomial vector coded\n");
						return true;
					}
					template <class GenericTruncator>
					void perform_hash_coded_multiplication(const GenericTruncator &trunc) {
						typedef coded_series_hash_table<cf_type1, max_fast_int> csht;
						typedef typename csht::term_type cterm;
						typedef typename csht::iterator c_iterator;
						csht cms;
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
						// TODO: rehash on m_retval here (since we know what the size is going to be)?
						// This would require the generic wrapper around the container of the series.
						term_type1 tmp_term;
						const c_iterator c_it_f = cms.end();
						for (c_iterator c_it = cms.begin(); c_it != c_it_f; ++c_it) {
							tmp_term.m_cf = c_it->m_cf;
							coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
							ancestor::m_retval.insert(tmp_term, ancestor::m_args_tuple);
						}
						__PDEBUG(std::cout << "Done polynomial hash coded\n");
					}
				private:
					// Temporary key used for the decodification in the truncator.
					// It is mutable because it is used as temporary decodification area.
					mutable typename term_type1::key_type	m_tmp_key;
			};
	};
}

#endif
