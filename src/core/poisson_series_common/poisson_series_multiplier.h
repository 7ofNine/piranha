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

#include <algorithm>
#include <boost/algorithm/minmax_element.hpp> // To calculate limits of multiplication.
#include <cstddef>
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_series_multiplier.h"
#include "../base_classes/null_truncator.h"
#include "../coded_series_hash_table.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../memory.h"
#include "../settings.h" // For debug.
#include "../type_traits.h" // For lightweight attribute.

namespace piranha
{
	/// Series multiplier specifically tuned for Poisson series.
	/**
	 * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
	// NOTE: isn't truncator's "accept" used here?
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
						this->find_input_min_max();
						calculate_result_min_max();
						this->determine_viability();
						// Use the selected truncator only if it really truncates, otherwise use the
						// null truncator.
						if (trunc.is_effective()) {
							ll_perform_multiplication(trunc);
						} else {
							if (is_lightweight<cf_type1>::value) {
								typedef typename ancestor::key_revlex_comparison key_revlex_comparison;
								std::sort(this->m_terms1.begin(),this->m_terms1.end(),key_revlex_comparison());
								std::sort(this->m_terms2.begin(),this->m_terms2.end(),key_revlex_comparison());
							}
							ll_perform_multiplication(null_truncator::template get_type<get_type>(*this));
						}
					}
				private:
					template <class GenericTruncator>
					void ll_perform_multiplication(const GenericTruncator &trunc) {
						if (!is_lightweight<cf_type1>::value || (this->m_terms1.size() < 10 && this->m_terms2.size() < 10)) {
							__PDEBUG(std::cout << "Heavy coefficient or small series, "
								"going for plain poisson series multiplication\n");
							this->perform_plain_multiplication(trunc);
						} else if (this->m_cr_is_viable) {
							this->code_keys();
							// We also need flavours here.
							cache_flavours();
							const term_type1 **t1 = &this->m_terms1[0];
							const term_type2 **t2 = &this->m_terms2[0];
							std::vector<cf_type1> cf1_cache;
							std::vector<cf_type2> cf2_cache;
							// NOTICE: this check is really compile-time, so we could probably avoid
							// having this "if" in favour of a meta-programmed chooser. However, the compiler
							// here probably just ditches this part while optimizing, so maybe it is really
							// not much worthwhile.
							if (is_lightweight<cf_type1>::value) {
								// Cache the values if cf_type1 is lightweight.
								const std::size_t size1 = this->m_size1, size2 = this->m_size2;
								cf1_cache.reserve(size1);
								cf2_cache.reserve(size2);
								for (std::size_t i = 0; i < size1; ++i) {
									cf1_cache.push_back(t1[i]->m_cf);
								}
								for (std::size_t i = 0; i < size2; ++i) {
									cf2_cache.push_back(t2[i]->m_cf);
								}
							}
							bool vec_res;
							if (this->is_sparse()) {
								vec_res = false;
							} else {
								if (is_lightweight<cf_type1>::value) {
									vec_res = perform_vector_coded_multiplication<ancestor::template cf_direct>(
										&cf1_cache[0],&cf2_cache[0],t1,t2,trunc);
								} else {
									vec_res = perform_vector_coded_multiplication<ancestor::template cf_from_term>(
										t1,t2,t1,t2,trunc);
								}
							}
							if (!vec_res) {
								__PDEBUG(std::cout << "Going for hash coded poisson series multiplication\n");
								if (is_lightweight<cf_type1>::value) {
									perform_hash_coded_multiplication<ancestor::template cf_direct>(
										&cf1_cache[0],&cf2_cache[0],t1,t2,trunc);
								} else {
									perform_hash_coded_multiplication<ancestor::template cf_from_term>(t1,t2,t1,t2,trunc);
								}
							}
						} else {
							__PDEBUG(std::cout << "Going for plain poisson series multiplication\n");
							this->perform_plain_multiplication(trunc);
						}
					}
					void calculate_result_min_max() {
						std::vector<mpz_class> tmp_vec(8);
						std::pair<typename std::vector<mpz_class>::const_iterator,
							std::vector<mpz_class>::const_iterator> min_max;
						const std::size_t size = this->m_size;
						for (std::size_t i = 0; i < size; ++i) {
							tmp_vec[0] = this->m_min_max1[i].second;
							tmp_vec[0] += this->m_min_max2[i].second;
							tmp_vec[1] = this->m_min_max1[i].first;
							tmp_vec[1] += this->m_min_max2[i].first;
							tmp_vec[2] = this->m_min_max1[i].second;
							tmp_vec[2] -= this->m_min_max2[i].first;
							tmp_vec[3] = this->m_min_max1[i].first;
							tmp_vec[3] -= this->m_min_max2[i].second;
							tmp_vec[4] = this->m_min_max1[i].first;
							tmp_vec[5] = this->m_min_max2[i].first;
							tmp_vec[6] = this->m_min_max1[i].second;
							tmp_vec[7] = this->m_min_max2[i].second;
							min_max = boost::minmax_element(tmp_vec.begin(), tmp_vec.end());
							this->m_res_min_max[i].first = *(min_max.first);
							this->m_res_min_max[i].second = *(min_max.second);
						}
						__PDEBUG(
							std::cout << "Mult limits are:\n";
						for (std::size_t i = 0; i < this->m_res_min_max.size(); ++i) {
						std::cout << this->m_res_min_max[i].first << ',' <<
							this->m_res_min_max[i].second << '\n';
						}
						);
					}
					// Store flavours of the series into own vectors.
					void cache_flavours() {
						std::size_t i;
						for (i = 0; i < this->m_size1; ++i) {
							m_flavours1[i] = this->m_terms1[i]->m_key.get_flavour();
						}
						for (i = 0; i < this->m_size2; ++i) {
							m_flavours2[i] = this->m_terms2[i]->m_key.get_flavour();
						}
					}
					struct vector_multiplier {
						vector_multiplier(const char *f1, const char *f2):m_f1(f1),m_f2(f2) {}
						template <template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
							class Term1, class Term2, class Ckey, class Trunc, class ResVec>
						bool run(
							const std::size_t &i, const std::size_t &j,
							const TermOrCf1 *tc1, const TermOrCf2 *tc2,
							const Term1 **t1, const Term2 **t2, const Ckey *ck1,
							const Ckey *ck2, const Trunc &trunc, ResVec *vc_res_pair, const ArgsTuple &args_tuple)
						{
							typedef CfGetter<cf_type1> get1;
							typedef CfGetter<cf_type2> get2;
							if (trunc.skip(&t1[i], &t2[j])) {
								return false;
							}
							// Cache values.
							cf_type1 *vc_res_cos = vc_res_pair->first, *vc_res_sin = vc_res_pair->second;
							const char *f1 = m_f1, *f2 = m_f2;
							// NOTE: Does it make sense here to define a method for coefficients like:
							// mult_by_and_insert_into<bool Sign>(cf2,retval,m_args_tuple)
							// so that we can avoid copying stuff around here and elsewhere?
							cf_type1 tmp_cf = get1::get(tc1[i]);
							tmp_cf.mult_by(get2::get(tc2[j]), args_tuple);
							const max_fast_int index_plus = ck1[i] + ck2[j], index_minus = ck1[i] - ck2[j];
							if (f1[i] == f2[j]) {
								if (f1[i]) {
									vc_res_cos[index_minus].add(tmp_cf, args_tuple);
									vc_res_cos[index_plus].add(tmp_cf, args_tuple);
								} else {
									vc_res_cos[index_minus].add(tmp_cf, args_tuple);
									vc_res_cos[index_plus].subtract(tmp_cf, args_tuple);
								}
							} else {
								if (f1[i]) {
									vc_res_sin[index_minus].subtract(tmp_cf, args_tuple);
									vc_res_sin[index_plus].add(tmp_cf, args_tuple);
								} else {
									vc_res_sin[index_minus].add(tmp_cf, args_tuple);
									vc_res_sin[index_plus].add(tmp_cf, args_tuple);
								}
							}
							return true;
						}
						const char	*m_f1;
						const char	*m_f2;
					};
					template <template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
						class Term1, class Term2, class GenericTruncator>
					bool perform_vector_coded_multiplication(const TermOrCf1 *tc1, const TermOrCf2 *tc2,
						const Term1 **t1, const Term2 **t2, const GenericTruncator &trunc)
					{
						std::vector<cf_type1,std_counting_allocator<cf_type1> > vc_cos, vc_sin;
						// Try to allocate the space for vector coded multiplication. We need two arrays of results,
						// one for cosines, one for sines.
						// The +1 is needed because we need the number of possible codes between min and max, e.g.:
						// coded_ancestor::m_h_min = 0, coded_ancestor::m_h_max = 2 --> n of codes = 3.
						piranha_assert(this->m_h_max - this->m_h_min + 1 >= 0);
						const std::size_t n_codes = static_cast<std::size_t>(this->m_h_max -
							this->m_h_min + 1);
						try {
							vc_cos.resize(n_codes);
							vc_sin.resize(n_codes);
						} catch (const std::bad_alloc &) {
							__PDEBUG(std::cout << "Not enough physical memory available for vector coded.\n");
							return false;
						} catch (const memory_error &) {
							__PDEBUG(std::cout << "Memory limit reached for vector coded.\n");
							return false;
						}
						__PDEBUG(std::cout << "Going for vector coded Poisson series multiplication\n");
						// Define the base pointers for storing the results of multiplication.
						// Please note that even if here it seems like we are going to write outside allocated memory,
						// the indices from the analysis of the coded series will prevent out-of-boundaries
						// reads/writes.
						const std::size_t size1 = this->m_size1, size2 = this->m_size2;
						const args_tuple_type &args_tuple = this->m_args_tuple;
						const max_fast_int *ck1 = &this->m_ckeys1[0], *ck2 = &this->m_ckeys2[0];
						std::pair<cf_type1 *, cf_type1 *> res(&vc_cos[0] - this->m_h_min, &vc_sin[0] - this->m_h_min);
						// Find out a suitable block size.
						static const std::size_t block_size =
							(2 << (ilg<isqrt<(settings::cache_size * 1024) / (sizeof(cf_type1))>::value>::value - 1));
						__PDEBUG(std::cout << "Block size: " << block_size << '\n');
						// Perform multiplication.
						vector_multiplier vm(&m_flavours1[0],&m_flavours2[0]);
						this->template blocked_multiplication<block_size,CfGetter>(
							size1,size2,tc1,tc2,t1,t2,ck1,ck2,trunc,&res,vm,args_tuple
						);
						__PDEBUG(std::cout << "Done multiplying\n");
						// Decode and insert the results into return value.
						cf_type1 *vc_res_cos = res.first, *vc_res_sin = res.second;
						term_type1 tmp_term;
						const max_fast_int i_f = this->m_h_max;
						for (max_fast_int i = this->m_h_min; i <= i_f; ++i) {
							vc_res_cos[i].divide_by(2,args_tuple);
							// Take a shortcut and check for ignorability of the coefficient here.
							// This way we avoid decodification, and all the series term insertion yadda-yadda.
							if (!vc_res_cos[i].is_ignorable(args_tuple)) {
								tmp_term.m_cf = vc_res_cos[i];
								this->decode(tmp_term.m_key, i);
								tmp_term.m_key.set_flavour(true);
								// Canonicalise in-place, so that we don't need to make further copies in the
								// main insertion function.
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						for (max_fast_int i = this->m_h_min; i <= i_f; ++i) {
							vc_res_sin[i].divide_by(2,args_tuple);
							if (!vc_res_sin[i].is_ignorable(args_tuple)) {
								tmp_term.m_cf = vc_res_sin[i];
								this->decode(tmp_term.m_key, i);
								tmp_term.m_key.set_flavour(false);
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done Poisson series vector coded\n");
						return true;
					}
					template <class Cterm>
					struct hash_multiplier {
						hash_multiplier(const char *f1, const char *f2):m_f1(f1),m_f2(f2) {}
						template <template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
							class Term1, class Term2, class Ckey, class Trunc, class HashSet>
						bool run(
							const std::size_t &i, const std::size_t &j,
							const TermOrCf1 *tc1, const TermOrCf2 *tc2,
							const Term1 **t1, const Term2 **t2, const Ckey *ck1,
							const Ckey *ck2, const Trunc &trunc, std::pair<HashSet *,HashSet *> *cms_pair,
							const ArgsTuple &args_tuple)
						{
							typedef CfGetter<cf_type1> get1;
							typedef CfGetter<cf_type2> get2;
							typedef typename HashSet::iterator c_iterator;
							if (trunc.skip(&t1[i], &t2[j])) {
								return false;
							}
							// Cache values.
							const char *f1 = m_f1, *f2 = m_f2;
							HashSet &cms_cos = *cms_pair->first, &cms_sin = *cms_pair->second;
							// TODO: here (and elsewhere, likely), we can avoid an extra copy by working with keys
							// and cfs instead of terms, generating only one coefficient and change its sign later
							// if needed - after insertion.
							// NOTE: cache tmp_term1 from external, as done in vector multiplier?
							Cterm tmp_term1(get1::get(tc1[i]), ck1[i]);
							// Handle the coefficient, with positive signs for now.
							tmp_term1.m_cf.mult_by(get2::get(tc2[j]), args_tuple);
							tmp_term1.m_ckey -= ck2[j];
							// Create the second term, using the first one's coefficient and the appropriate code.
							Cterm tmp_term2(tmp_term1.m_cf, ck1[i] + ck2[j]);
							// Now fix flavours and coefficient signs.
							if (f1[i] == f2[j]) {
								if (!f1[i]) {
									tmp_term2.m_cf.invert_sign(args_tuple);
								}
								// Insert into cosine container.
								std::pair<bool,c_iterator> res = cms_cos.find(tmp_term1);
								if (res.first) {
									res.second->m_cf.add(tmp_term1.m_cf, args_tuple);
								} else {
									cms_cos.insert(tmp_term1,res.second);
								}
								res = cms_cos.find(tmp_term2);
								if (res.first) {
									res.second->m_cf.add(tmp_term2.m_cf, args_tuple);
								} else {
									cms_cos.insert(tmp_term2,res.second);
								}
							} else {
								if (f1[i]) {
									tmp_term1.m_cf.invert_sign(args_tuple);
								}
								// Insert into sine container.
								std::pair<bool,c_iterator> res = cms_sin.find(tmp_term1);
								if (res.first) {
									res.second->m_cf.add(tmp_term1.m_cf, args_tuple);
								} else {
									cms_sin.insert(tmp_term1,res.second);
								}
								res = cms_sin.find(tmp_term2);
								if (res.first) {
									res.second->m_cf.add(tmp_term2.m_cf, args_tuple);
								} else {
									cms_sin.insert(tmp_term2,res.second);
								}
							}
							return true;
						}
						const char	*m_f1;
						const char	*m_f2;
					};
					template <template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
						class Term1, class Term2, class GenericTruncator>
					void perform_hash_coded_multiplication(const TermOrCf1 *tc1, const TermOrCf2 *tc2,
						const Term1 **t1, const Term2 **t2, const GenericTruncator &trunc)
					{
						typedef typename coded_ancestor::template coded_term_type<cf_type1,max_fast_int> cterm;
						typedef coded_series_hash_table<cterm, std::allocator<char> > csht;
						typedef typename csht::iterator c_iterator;
						// Let's find a sensible size hint.
						const std::size_t size_hint = static_cast<std::size_t>(
							std::max<double>(this->m_density1,this->m_density2) * this->m_h_tot);
						csht cms_cos(size_hint), cms_sin(size_hint);
						std::pair<csht *, csht *> res(&cms_cos,&cms_sin);
						const std::size_t size1 = this->m_size1, size2 = this->m_size2;
						const args_tuple_type &args_tuple = this->m_args_tuple;
						const max_fast_int *ck1 = &this->m_ckeys1[0], *ck2 = &this->m_ckeys2[0];
						// Find out a suitable block size.
						static const std::size_t block_size =
							(2 << (ilg<isqrt<(settings::cache_size * 1024) / (sizeof(cterm))>::value>::value - 1));
						__PDEBUG(std::cout << "Block size: " << block_size << '\n');
						hash_multiplier<cterm> hm(&m_flavours1[0],&m_flavours2[0]);
						this->template blocked_multiplication<block_size,CfGetter>(
							size1,size2,tc1,tc2,t1,t2,ck1,ck2,trunc,&res,hm,args_tuple
						);
						__PDEBUG(std::cout << "Done Poisson series hash coded multiplying\n");
						term_type1 tmp_term;
						{
							const c_iterator c_it_f = cms_cos.end();
							for (c_iterator c_it = cms_cos.begin(); c_it != c_it_f; ++c_it) {
								tmp_term.m_cf = c_it->m_cf;
								tmp_term.m_cf.divide_by(2,args_tuple);
								this->decode(tmp_term.m_key, c_it->m_ckey);
								tmp_term.m_key.set_flavour(true);
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						{
							const c_iterator c_it_f = cms_sin.end();
							for (c_iterator c_it = cms_sin.begin(); c_it != c_it_f; ++c_it) {
								tmp_term.m_cf = c_it->m_cf;
								tmp_term.m_cf.divide_by(2,args_tuple);
								this->decode(tmp_term.m_key, c_it->m_ckey);
								tmp_term.m_key.set_flavour(false);
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done Poisson series hash coded\n");
					}
				private:
					// For Poisson series we also need flavours.
					std::vector<char>	m_flavours1;
					std::vector<char>	m_flavours2;
			};
	};
}

#endif
