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
#include <boost/thread/thread.hpp>
#include <cmath>
#include <cstddef>
#include <exception>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_series_multiplier.h"
#include "../base_classes/null_truncator.h"
#include "../coded_series_hash_table.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../mp.h"
#include "../runtime.h"
#include "../settings.h" // For debug and cache size.
#include "../type_traits.h" // For lightweight attribute.

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
				public:
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<Series1,Series2,ArgsTuple> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
						ancestor(s1, s2, retval, args_tuple) {}
					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication()
					{
// 						const std::size_t n = this->m_ranges.size();
// 						piranha_assert(n > 0);
// 						if (n == 1) {
// 							worker(this,&(this->m_retval),0);
// 							worker();
// 						} else {
// 							boost::thread_group() tg;
// 							std::vector<Series1> retvals(n,Series1());
// 							for (std::size_t i = 0; i < n; ++i) {
// 								tg.create_thread(worker(this,&retvals[i],i));
// 							}
// 							tg.join_all();
// 							// Take the retvals and insert them into final retval.
// 							// TODO
// 						}
						// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
						// coefficients and coded keys.
						const truncator_type trunc(this->m_terms1,this->m_terms2,this->m_args_tuple);
						this->find_input_min_max();
						calculate_result_min_max();
						this->determine_viability();
						if (trunc.is_effective()) {
							ll_perform_multiplication(trunc);
						} else {
							// NOTE: maybe here it is worth to sort them anyway, given the fact
							// that ascending codes should still be beneficial for cache usage.
							// We want to sort them this way if we are not truncating and
							// coefficients are lightweight, in order to optimize cache memory usage
							if (is_lightweight<cf_type1>::value) {
								typedef typename ancestor::key_revlex_comparison key_revlex_comparison;
								std::sort(this->m_terms1.begin(),this->m_terms1.end(),key_revlex_comparison());
								std::sort(this->m_terms2.begin(),this->m_terms2.end(),key_revlex_comparison());
							}
							ll_perform_multiplication(null_truncator::template get_type<Series1,Series2,ArgsTuple>(
								this->m_terms1,this->m_terms2,this->m_args_tuple
							));
						}
					}
				private:
					template <class GenericTruncator>
					void ll_perform_multiplication(const GenericTruncator &trunc) {
						// TODO: these checks maybe should go earlier in order to avoid redoing truncator check in plain multiplication?
						if (!is_lightweight<cf_type1>::value || (this->m_terms1.size() < 10 && this->m_terms2.size() < 10)) {
							__PDEBUG(std::cout << "Heavy coefficient or small series, "
								"going for plain polynomial multiplication\n");
							this->perform_plain_multiplication();
						} else if (this->m_cr_is_viable) {
							this->code_keys();
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
								__PDEBUG(std::cout << "Going for hash coded polynomial multiplication\n");
								if (is_lightweight<cf_type1>::value) {
									perform_hash_coded_multiplication<ancestor::template cf_direct>(
										&cf1_cache[0],&cf2_cache[0],t1,t2,trunc);
								} else {
									perform_hash_coded_multiplication<ancestor::template cf_from_term>(t1,t2,t1,t2,trunc);
								}
							}
						} else {
							__PDEBUG(std::cout << "Going for plain polynomial multiplication\n");
							this->perform_plain_multiplication();
						}
					}
					// TODO: better rename result_min_max here and everywhere. It is not really the min/max
					// values for the result, because the min/max values of the input series are also taken
					// into account to establish the codification parameters for the input series _and_
					// the resulting series.
					void calculate_result_min_max() {
						std::vector<mp_integer> tmp_vec(6);
						std::pair<typename std::vector<mp_integer>::const_iterator,
							std::vector<mp_integer>::const_iterator> min_max;
						for (std::size_t i = 0; i < this->m_size; ++i) {
							tmp_vec[0] = this->m_min_max1[i].second;
							tmp_vec[0] += this->m_min_max2[i].second;
							tmp_vec[1] = this->m_min_max1[i].first;
							tmp_vec[1] += this->m_min_max2[i].first;
							tmp_vec[2] = this->m_min_max1[i].first;
							tmp_vec[3] = this->m_min_max2[i].first;
							tmp_vec[4] = this->m_min_max1[i].second;
							tmp_vec[5] = this->m_min_max2[i].second;
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
					template <template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
						class Ckey, class Trunc>
					struct vector_functor {
						vector_functor(const TermOrCf1 *tc1, const TermOrCf2 *tc2,
							const term_type1 **t1, const term_type2 **t2,
							const Ckey *ck1, const Ckey *ck2,
							const Trunc &trunc, cf_type1 *vc_res, const ArgsTuple &args_tuple):
							m_tc1(tc1),m_tc2(tc2),m_t1(t1),m_t2(t2),m_ck1(ck1),m_ck2(ck2),m_trunc(trunc),
							m_vc_res(vc_res),m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							typedef CfGetter<cf_type1> get1;
							typedef CfGetter<cf_type2> get2;
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							// Calculate index of the result.
							const max_fast_int res_index = m_ck1[i] + m_ck2[j];
							if (m_trunc.accept(res_index)) {
								m_vc_res[res_index].addmul(get1::get(m_tc1[i]), get2::get(m_tc2[j]), m_args_tuple);
							}
							return true;
						}
						const TermOrCf1		*m_tc1;
						const TermOrCf2		*m_tc2;
						const term_type1	**m_t1;
						const term_type2	**m_t2;
						const Ckey		*m_ck1;
						const Ckey		*m_ck2;
						const Trunc		&m_trunc;
						cf_type1		*m_vc_res;
						const ArgsTuple		&m_args_tuple;
					};
					template <template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
						class Term1, class Term2, class GenericTruncator>
					bool perform_vector_coded_multiplication(const TermOrCf1 *tc1, const TermOrCf2 *tc2,
						const Term1 **t1, const Term2 **t2, const GenericTruncator &trunc)
					{
						std::vector<cf_type1,std_counting_allocator<cf_type1> > vc;
						// Try to allocate the space for vector coded multiplication.
						// The +1 is needed because we need the number of possible codes between min and max, e.g.:
						// coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
						piranha_assert(this->m_h_max - this->m_h_min + 1 >= 0);
						const std::size_t n_codes = static_cast<std::size_t>(this->m_h_max - this->m_h_min + 1);
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
						// Please note that even if here it seems like we are going to write outside allocated memory,
						// the indices from the analysis of the coded series will prevent out-of-boundaries
						// reads/writes.
						const std::size_t size1 = this->m_size1, size2 = this->m_size2;
						const max_fast_int *ck1 = &this->m_ckeys1[0], *ck2 = &this->m_ckeys2[0];
						const args_tuple_type &args_tuple = this->m_args_tuple;
						cf_type1 *vc_res =  &vc[0] - this->m_h_min;
						// Find out a suitable block size.
						const std::size_t block_size = 2 <<
							((std::size_t)log2(std::max(16.,std::sqrt((settings::cache_size * 1024) / (sizeof(cf_type1) * runtime::get_n_cur_threads())))) - 1);
						__PDEBUG(std::cout << "Block size: " << block_size << '\n');
						// Perform multiplication.
						vector_functor<CfGetter,TermOrCf1,TermOrCf2,max_fast_int,GenericTruncator>
							vm(tc1,tc2,t1,t2,ck1,ck2,trunc,vc_res,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,vm);
						__PDEBUG(std::cout << "Done multiplying\n");
						const max_fast_int i_f = this->m_h_max;
						// Decode and insert the results into return value.
						term_type1 tmp_term;
						for (max_fast_int i = this->m_h_min; i <= i_f; ++i) {
							// Take a shortcut and check for ignorability of the coefficient here.
							// This way we avoid decodification, and all the series term insertion yadda-yadda.
							// NOTE: wouldn't it be better if insert() were smart enough to do these checks first
							// and reduce its workload?
							if (!vc_res[i].is_ignorable(args_tuple)) {
								tmp_term.m_cf = vc_res[i];
								coded_ancestor::decode(tmp_term.m_key, i);
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						__PDEBUG(std::cout << "Done polynomial vector coded.\n");
						return true;
					}
					template <class Cterm,template <class> class CfGetter, class TermOrCf1, class TermOrCf2,
						class Ckey, class Trunc, class HashSet>
					struct hash_functor {
						hash_functor(Cterm &cterm, const TermOrCf1 *tc1, const TermOrCf2 *tc2,
							const term_type1 **t1, const term_type2 **t2,
							const Ckey *ck1, const Ckey *ck2,
							const Trunc &trunc, HashSet *cms, const ArgsTuple &args_tuple):
							m_cterm(cterm),m_tc1(tc1),m_tc2(tc2),m_t1(t1),m_t2(t2),m_ck1(ck1),m_ck2(ck2),
							m_trunc(trunc),m_cms(cms),m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							typedef CfGetter<cf_type1> get1;
							typedef CfGetter<cf_type2> get2;
							typedef typename HashSet::iterator c_iterator;
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							m_cterm.m_ckey = m_ck1[i];
							m_cterm.m_ckey += m_ck2[j];
							if (m_trunc.accept(m_cterm.m_ckey)) {
								const std::pair<bool,c_iterator> res = m_cms->find(m_cterm);
								if (res.first) {
									res.second->m_cf.addmul(get1::get(m_tc1[i]), get2::get(m_tc2[j]), m_args_tuple);
								} else {
									// Assign to the temporary term the old cf (new_key is already assigned).
									m_cterm.m_cf = get1::get(m_tc1[i]);
									// Multiply the old term by the second term.
									m_cterm.m_cf.mult_by(get2::get(m_tc2[j]), m_args_tuple);
									m_cms->insert(m_cterm,res.second);
								}
							}
							return true;
						}
						Cterm			&m_cterm;
						const TermOrCf1		*m_tc1;
						const TermOrCf2		*m_tc2;
						const term_type1	**m_t1;
						const term_type2	**m_t2;
						const Ckey		*m_ck1;
						const Ckey		*m_ck2;
						const Trunc		&m_trunc;
						HashSet			*m_cms;
						const ArgsTuple		&m_args_tuple;
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
						const std::size_t size1 = this->m_size1, size2 = this->m_size2;
						const max_fast_int *ck1 = &this->m_ckeys1[0], *ck2 = &this->m_ckeys2[0];
						const args_tuple_type &args_tuple = this->m_args_tuple;
						csht cms(size_hint);
						// Find out a suitable block size.
						const std::size_t block_size = 2 <<
							((std::size_t)log2(std::max(16.,std::sqrt((settings::cache_size * 1024) / (sizeof(cterm) * runtime::get_n_cur_threads())))) - 1);
						__PDEBUG(std::cout << "Block size: " << block_size << '\n');
						cterm tmp_cterm;
						hash_functor<cterm,CfGetter,TermOrCf1,TermOrCf2,max_fast_int,GenericTruncator,csht>
							hm(tmp_cterm,tc1,tc2,t1,t2,ck1,ck2,trunc,&cms,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,hm);
						__PDEBUG(std::cout << "Done polynomial hash coded multiplying\n");
						// Decode and insert into retval.
						term_type1 tmp_term;
						// TODO: add debug info about cms' size here.
						const c_iterator c_it_f = cms.end();
						for (c_iterator c_it = cms.begin(); c_it != c_it_f; ++c_it) {
							tmp_term.m_cf = c_it->m_cf;
							coded_ancestor::decode(tmp_term.m_key, c_it->m_ckey);
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
