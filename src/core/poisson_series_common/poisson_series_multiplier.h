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
#include <cmath>
#include <cstddef>
#include <exception>
#include <iterator>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_series_multiplier.h"
#include "../base_classes/null_truncator.h"
#include "../coded_series_hash_table.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"
#include "../memory.h"
#include "../settings.h" // For debug.
#include "../type_traits.h" // For lightweight attribute.
#include "../utils.h"

namespace piranha
{
	/// Series multiplier specifically tuned for Poisson series.
	/**
	 * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
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
					typedef typename Truncator::template get_type<Series1,Series2,ArgsTuple> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
						ancestor(s1, s2, retval, args_tuple),
						m_flavours1(ancestor::m_size1), m_flavours2(ancestor::m_size2) {}
					/// Perform multiplication and place the result into m_retval.
					void perform_multiplication()
					{
						// Cache term pointers.
						utils::cache_terms_pointers(this->m_s1,this->m_terms1);
						utils::cache_terms_pointers(this->m_s2,this->m_terms2);
						// NOTE: hard coded value of 1000.
						// NOTE: share in a coded multiplier toolbox?
						if (!is_lightweight<cf_type1>::value || double(this->m_size1) * double(this->m_size2) < 1000) {
							__PDEBUG(std::cout << "Heavy coefficient or small series, "
								"going for plain multiplication\n");
							this->perform_plain_multiplication();
							return;
						}
						// Build the truncator here, _before_ coding. Otherwise we mess up the relation between
						// coefficients and coded keys.
						const truncator_type trunc(this->m_terms1,this->m_terms2,this->m_args_tuple);
						this->find_input_min_max();
						calculate_result_min_max();
						this->determine_viability();
						if (!this->m_cr_is_viable) {
							__PDEBUG(std::cout << "Series not suitable for coded representation, going for plain multiplication\n");
							this->perform_plain_multiplication();
							return;
						}
						if (trunc.is_effective()) {
							ll_perform_multiplication(trunc);
						} else {
							// Sort input series for better cache usage and multi-threaded implementation.
							typedef typename ancestor::key_revlex_comparison key_revlex_comparison;
							std::sort(this->m_terms1.begin(),this->m_terms1.end(),key_revlex_comparison());
							std::sort(this->m_terms2.begin(),this->m_terms2.end(),key_revlex_comparison());
							ll_perform_multiplication(null_truncator::template get_type<Series1,Series2,ArgsTuple>(
								this->m_terms1,this->m_terms2,this->m_args_tuple
							));
						}
					}
				private:
					template <class GenericTruncator>
					void ll_perform_multiplication(const GenericTruncator &trunc)
					{
						this->code_keys();
						// We also need flavours here.
						cache_flavours();
						const term_type1 **t1 = &this->m_terms1[0];
						const term_type2 **t2 = &this->m_terms2[0];
						// Cache the coefficients.
						// NOTE: c++0x lambdas here.
						std::vector<cf_type1> cf1_cache;
						std::vector<cf_type2> cf2_cache;
						std::insert_iterator<std::vector<cf_type1> > i_it1(cf1_cache,cf1_cache.begin());
						std::insert_iterator<std::vector<cf_type2> > i_it2(cf2_cache,cf2_cache.begin());
						std::transform(this->m_terms1.begin(),this->m_terms1.end(),i_it1,typename ancestor::template ptr_cf_extractor<term_type1>());
						std::transform(this->m_terms2.begin(),this->m_terms2.end(),i_it2,typename ancestor::template ptr_cf_extractor<term_type2>());
						bool vec_res;
						if (this->is_sparse()) {
							vec_res = false;
						} else {
							vec_res = perform_vector_coded_multiplication(&cf1_cache[0],&cf2_cache[0],t1,t2,trunc);
						}
						if (!vec_res) {
							__PDEBUG(std::cout << "Going for hash coded poisson series multiplication\n");
							perform_hash_coded_multiplication(&cf1_cache[0],&cf2_cache[0],t1,t2,trunc);
						}
					}
					void calculate_result_min_max()
					{
						std::vector<mp_integer> tmp_vec(8);
						std::pair<typename std::vector<mp_integer>::const_iterator,
							std::vector<mp_integer>::const_iterator> min_max;
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
					void cache_flavours()
					{
						std::size_t i;
						for (i = 0; i < this->m_size1; ++i) {
							m_flavours1[i] = this->m_terms1[i]->m_key.get_flavour();
						}
						for (i = 0; i < this->m_size2; ++i) {
							m_flavours2[i] = this->m_terms2[i]->m_key.get_flavour();
						}
					}
					template <class GenericTruncator>
					struct vector_functor {
						vector_functor(const char *f1, const char *f2,
							const cf_type1 *tc1, const cf_type2 *tc2,
							const term_type1 **t1, const term_type2 **t2,
							const max_fast_int *ck1, const max_fast_int *ck2,
							const GenericTruncator &trunc, std::pair<cf_type1 *, cf_type1 *> *vc_res_pair, const ArgsTuple &args_tuple):
							m_f1(f1),m_f2(f2),m_tc1(tc1),m_tc2(tc2),m_t1(t1),m_t2(t2),m_ck1(ck1),m_ck2(ck2),m_trunc(trunc),
							m_vc_res_pair(vc_res_pair),m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							// Cache values.
							cf_type1 *vc_res_cos = m_vc_res_pair->first, *vc_res_sin = m_vc_res_pair->second;
							const char *f1 = m_f1, *f2 = m_f2;
							// NOTE: Does it make sense here to define a method for coefficients like:
							// mult_by_and_insert_into<bool Sign>(cf2,retval,m_args_tuple)
							// so that we can avoid copying stuff around here and elsewhere?
							cf_type1 tmp_cf = m_tc1[i];
							tmp_cf.mult_by(m_tc2[j], m_args_tuple);
							const max_fast_int index_plus = m_ck1[i] + m_ck2[j], index_minus = m_ck1[i] - m_ck2[j];
							if (f1[i] == f2[j]) {
								if (f1[i]) {
									vc_res_cos[index_minus].add(tmp_cf, m_args_tuple);
									vc_res_cos[index_plus].add(tmp_cf, m_args_tuple);
								} else {
									vc_res_cos[index_minus].add(tmp_cf, m_args_tuple);
									vc_res_cos[index_plus].subtract(tmp_cf, m_args_tuple);
								}
							} else {
								if (f1[i]) {
									vc_res_sin[index_minus].subtract(tmp_cf, m_args_tuple);
									vc_res_sin[index_plus].add(tmp_cf, m_args_tuple);
								} else {
									vc_res_sin[index_minus].add(tmp_cf, m_args_tuple);
									vc_res_sin[index_plus].add(tmp_cf, m_args_tuple);
								}
							}
							return true;
						}
						const char				*m_f1;
						const char				*m_f2;
						const cf_type1				*m_tc1;
						const cf_type2				*m_tc2;
						const term_type1			**m_t1;
						const term_type2			**m_t2;
						const max_fast_int			*m_ck1;
						const max_fast_int			*m_ck2;
						const GenericTruncator			&m_trunc;
						std::pair<cf_type1 *, cf_type1 *>	*m_vc_res_pair;
						const ArgsTuple				&m_args_tuple;
					};
					template <class GenericTruncator>
					bool perform_vector_coded_multiplication(const cf_type1 *tc1, const cf_type2 *tc2,
						const term_type1 **t1, const term_type2 **t2, const GenericTruncator &trunc)
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
						const std::size_t block_size = 2 <<
							((std::size_t)log2(std::max(16.,std::sqrt((settings::cache_size * 1024) / (sizeof(cf_type1))))) - 1);
						__PDEBUG(std::cout << "Block size: " << block_size << '\n';)
						// Perform multiplication.
						vector_functor<GenericTruncator> vm(&m_flavours1[0],&m_flavours2[0],tc1,tc2,t1,t2,ck1,ck2,trunc,&res,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,vm);
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
					template <class Cterm, class Ckey, class GenericTruncator, class HashSet>
					struct hash_functor {
						hash_functor(const char *f1, const char *f2,
							const cf_type1 *tc1, const cf_type2 *tc2,
							const term_type1 **t1, const term_type2 **t2,
							const Ckey *ck1, const Ckey *ck2,
							const GenericTruncator &trunc, std::pair<HashSet *,HashSet *> *cms, const ArgsTuple &args_tuple):
							m_f1(f1),m_f2(f2),
							m_tc1(tc1),m_tc2(tc2),m_t1(t1),m_t2(t2),m_ck1(ck1),m_ck2(ck2),
							m_trunc(trunc),m_cms(cms),m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							typedef typename HashSet::iterator c_iterator;
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							// Cache values.
							const char *f1 = m_f1, *f2 = m_f2;
							HashSet &cms_cos = *m_cms->first, &cms_sin = *m_cms->second;
							// TODO: here (and elsewhere, likely), we can avoid an extra copy by working with keys
							// and cfs instead of terms, generating only one coefficient and change its sign later
							// if needed - after insertion.
							// NOTE: cache tmp_term1 from external, as done in vector multiplier?
							Cterm tmp_term1(m_tc1[i], m_ck1[i]);
							// Handle the coefficient, with positive signs for now.
							tmp_term1.m_cf.mult_by(m_tc2[j], m_args_tuple);
							tmp_term1.m_ckey -= m_ck2[j];
							// Create the second term, using the first one's coefficient and the appropriate code.
							Cterm tmp_term2(tmp_term1.m_cf, m_ck1[i] + m_ck2[j]);
							// Now fix flavours and coefficient signs.
							if (f1[i] == f2[j]) {
								if (!f1[i]) {
									tmp_term2.m_cf.invert_sign(m_args_tuple);
								}
								// Insert into cosine container.
								std::pair<bool,c_iterator> res = cms_cos.find(tmp_term1);
								if (res.first) {
									res.second->m_cf.add(tmp_term1.m_cf, m_args_tuple);
								} else {
									cms_cos.insert(tmp_term1,res.second);
								}
								res = cms_cos.find(tmp_term2);
								if (res.first) {
									res.second->m_cf.add(tmp_term2.m_cf, m_args_tuple);
								} else {
									cms_cos.insert(tmp_term2,res.second);
								}
							} else {
								if (f1[i]) {
									tmp_term1.m_cf.invert_sign(m_args_tuple);
								}
								// Insert into sine container.
								std::pair<bool,c_iterator> res = cms_sin.find(tmp_term1);
								if (res.first) {
									res.second->m_cf.add(tmp_term1.m_cf, m_args_tuple);
								} else {
									cms_sin.insert(tmp_term1,res.second);
								}
								res = cms_sin.find(tmp_term2);
								if (res.first) {
									res.second->m_cf.add(tmp_term2.m_cf, m_args_tuple);
								} else {
									cms_sin.insert(tmp_term2,res.second);
								}
							}
							return true;
						}
						const char			*m_f1;
						const char			*m_f2;
						const cf_type1			*m_tc1;
						const cf_type2			*m_tc2;
						const term_type1		**m_t1;
						const term_type2		**m_t2;
						const Ckey			*m_ck1;
						const Ckey			*m_ck2;
						const GenericTruncator		&m_trunc;
						std::pair<HashSet *,HashSet *>	*m_cms;
						const ArgsTuple			&m_args_tuple;
					};
					template <class GenericTruncator>
					void perform_hash_coded_multiplication(const cf_type1 *tc1, const cf_type2 *tc2,
						const term_type1 **t1, const term_type2 **t2, const GenericTruncator &trunc)
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
						const std::size_t block_size = 2 <<
							((std::size_t)log2(std::max(16.,std::sqrt((settings::cache_size * 1024) / (sizeof(cterm))))) - 1);
						__PDEBUG(std::cout << "Block size: " << block_size << '\n';)
						hash_functor<cterm,max_fast_int,GenericTruncator,csht>
							hm(&m_flavours1[0],&m_flavours2[0],tc1,tc2,t1,t2,ck1,ck2,trunc,&res,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,hm);
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
