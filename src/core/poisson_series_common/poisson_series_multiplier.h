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

#include <boost/lambda/lambda.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <cstddef>
#include <exception>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/base_series_multiplier.h"
#include "../base_classes/coded_multiplier.h"
#include "../coded_hash_table.h"
#include "../config.h" // For p_static_check.
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../settings.h" // For debug.
#include "../stats.h"

namespace piranha
{
	// Select types of operations on coded series codes depending on whether we are dealing with a Poisson or Fourier series.
	template <int EchelonLevel>
	struct poisson_series_multiplier_ops_selector
	{
		p_static_check(EchelonLevel == 0,"");
		typedef boost::tuple<boost::false_type> type;
	};

	template <>
	struct poisson_series_multiplier_ops_selector<1>
	{
		typedef boost::tuple<boost::false_type,boost::true_type> type;
	};

	/// Series multiplier specifically tuned for Poisson series.
	/**
	 * This multiplier internally will use coded arithmetics if possible, otherwise it will operate just
	 * like piranha::base_series_multiplier.
	 */
	class poisson_series_multiplier
	{
		public:
			template <class Series1, class Series2, class ArgsTuple, class Truncator>
			class get_type:
				public base_series_multiplier<Series1, Series2, ArgsTuple, Truncator,
				get_type<Series1, Series2, ArgsTuple, Truncator> > ,
				public coded_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,
				typename poisson_series_multiplier_ops_selector<Series1::echelon_level>::type>
			{
					typedef base_series_multiplier< Series1, Series2, ArgsTuple, Truncator,
						get_type<Series1, Series2, ArgsTuple, Truncator> > ancestor;
					typedef coded_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,
						typename poisson_series_multiplier_ops_selector<Series1::echelon_level>::type> coded_ancestor;
					friend class coded_multiplier<get_type<Series1, Series2, ArgsTuple, Truncator>,Series1,Series2,
						typename poisson_series_multiplier_ops_selector<Series1::echelon_level>::type>;
					typedef typename ancestor::term_type1 term_type1;
					typedef typename ancestor::term_type2 term_type2;
					typedef typename final_cf<Series1>::type cf_type1;
					typedef typename final_cf<Series2>::type cf_type2;
				public:
					typedef Series1 series_type1;
					typedef Series2 series_type2;
					typedef ArgsTuple args_tuple_type;
					typedef typename Truncator::template get_type<Series1,Series2,ArgsTuple> truncator_type;
					get_type(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
						ancestor(s1, s2, retval, args_tuple) {}
					template <class GenericTruncator>
					void ll_perform_multiplication(const GenericTruncator &trunc)
					{
						// We also need flavours here.
						cache_flavours();
						// Procede with method from coded ancestor.
						coded_ancestor::ll_perform_multiplication(trunc);
					}
					// Store flavours of the series into own vectors.
					void cache_flavours()
					{
						m_flavours1.resize(this->m_terms1.size());
						m_flavours2.resize(this->m_terms2.size());
						std::size_t i;
						for (i = 0; i < this->m_terms1.size(); ++i) {
							m_flavours1[i] = this->m_terms1[i]->m_key.get_flavour();
						}
						for (i = 0; i < this->m_terms2.size(); ++i) {
							m_flavours2[i] = this->m_terms2[i]->m_key.get_flavour();
						}
					}
					template <class GenericTruncator>
					struct vector_functor {
						vector_functor(const char *f1, const char *f2,
							const cf_type1 *tc1, const cf_type2 *tc2,
							const term_type1 **t1, const term_type2 **t2,
							const max_fast_int *ck1, const max_fast_int *ck2a, const max_fast_int *ck2b,
							const GenericTruncator &trunc, std::pair<cf_type1 *, cf_type1 *> *vc_res_pair, const ArgsTuple &args_tuple):
							m_f1(f1),m_f2(f2),m_tc1(tc1),m_tc2(tc2),m_t1(t1),m_t2(t2),m_ck1(ck1),m_ck2a(ck2a),m_ck2b(ck2b),m_trunc(trunc),
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
							const max_fast_int index_plus = m_ck1[i] + m_ck2a[j], index_minus = m_ck1[i] + m_ck2b[j];
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
						const max_fast_int			*m_ck2a;
						const max_fast_int			*m_ck2b;
						const GenericTruncator			&m_trunc;
						std::pair<cf_type1 *, cf_type1 *>	*m_vc_res_pair;
						const ArgsTuple				&m_args_tuple;
					};
					template <class GenericTruncator>
					bool perform_vector_coded_multiplication(const cf_type1 *tc1, const cf_type2 *tc2,
						const term_type1 **t1, const term_type2 **t2, const GenericTruncator &trunc)
					{
						stats::trace_stat("mult_st",std::size_t(0),boost::lambda::_1 + 1);
						std::vector<cf_type1,std_counting_allocator<cf_type1> > vc_cos, vc_sin;
						// Try to allocate the space for vector coded multiplication. We need two arrays of results,
						// one for cosines, one for sines.
						// The +1 is needed because we need the number of possible codes between min and max, e.g.:
						// coded_ancestor::m_h_min = 0, coded_ancestor::m_h_max = 2 --> n of codes = 3.
						piranha_assert(boost::numeric::width(this->m_fast_h) + 1 >= 0);
						const std::size_t n_codes = boost::numeric_cast<std::size_t>(boost::numeric::width(this->m_fast_h) + 1);
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
						const std::size_t size1 = this->m_terms1.size(), size2 = this->m_terms2.size();
						piranha_assert(size1 && size2);
						const max_fast_int *ck1 = &this->m_ckeys1[0], *ck2a = &this->m_ckeys2a[0], *ck2b = &this->m_ckeys2b[0];
						const args_tuple_type &args_tuple = this->m_args_tuple;
						std::pair<cf_type1 *, cf_type1 *> res(&vc_cos[0] - this->m_fast_h.lower(), &vc_sin[0] - this->m_fast_h.lower());
						// Find out a suitable block size.
						const std::size_t block_size = this->template compute_block_size<sizeof(cf_type1)>();
						__PDEBUG(std::cout << "Block size: " << block_size << '\n';)
						// Perform multiplication.
						vector_functor<GenericTruncator> vm(&m_flavours1[0],&m_flavours2[0],tc1,tc2,t1,t2,ck1,ck2a,ck2b,trunc,&res,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,vm);
						__PDEBUG(std::cout << "Done multiplying\n");
						// Decode and insert the results into return value.
						cf_type1 *vc_res_cos = res.first, *vc_res_sin = res.second;
						term_type1 tmp_term;
						const max_fast_int i_f = this->m_fast_h.upper();
						for (max_fast_int i = this->m_fast_h.lower(); i <= i_f; ++i) {
							vc_res_cos[i].divide_by(2,args_tuple);
							// Take a shortcut and check for ignorability of the coefficient here.
							// This way we avoid decodification, and all the series term insertion yadda-yadda.
							if (!vc_res_cos[i].is_ignorable(args_tuple)) {
								this->decode(vc_res_cos[i],i,tmp_term);
								tmp_term.m_key.set_flavour(true);
								// Canonicalise in-place, so that we don't need to make further copies in the
								// main insertion function.
								if (!tmp_term.is_canonical(args_tuple)) {
									tmp_term.canonicalise(args_tuple);
								}
								this->m_retval.insert(tmp_term, args_tuple);
							}
						}
						for (max_fast_int i = this->m_fast_h.lower(); i <= i_f; ++i) {
							vc_res_sin[i].divide_by(2,args_tuple);
							if (!vc_res_sin[i].is_ignorable(args_tuple)) {
								this->decode(vc_res_sin[i],i,tmp_term);
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
							const Ckey *ck1, const Ckey *ck2a, const Ckey *ck2b,
							const GenericTruncator &trunc, std::pair<HashSet *,HashSet *> *cms, Cterm *tmp_term1, Cterm *tmp_term2,
							const ArgsTuple &args_tuple):
							m_f1(f1),m_f2(f2),
							m_tc1(tc1),m_tc2(tc2),m_t1(t1),m_t2(t2),m_ck1(ck1),m_ck2a(ck2a),m_ck2b(ck2b),
							m_trunc(trunc),m_cms(cms),
							m_tmp_term1(tmp_term1),m_tmp_term2(tmp_term2),
							m_args_tuple(args_tuple) {}
						bool operator()(const std::size_t &i, const std::size_t &j)
						{
							typedef typename HashSet::iterator c_iterator;
							if (m_trunc.skip(&m_t1[i], &m_t2[j])) {
								return false;
							}
							// Cache values.
							const char *f1 = m_f1, *f2 = m_f2;
							HashSet &cms_cos = *m_cms->first, &cms_sin = *m_cms->second;
							// NOTE: here (and elsewhere, likely), we can avoid an extra copy by working with keys
							// and cfs instead of terms, generating only one coefficient and change its sign later
							// if needed - after insertion <-- not sure this comment is still relevant....
							Cterm &tmp_term1 = *m_tmp_term1;
							tmp_term1.first = m_tc1[i];
							tmp_term1.second = m_ck1[i];
							// Handle the coefficient, with positive signs for now.
							tmp_term1.first.mult_by(m_tc2[j], m_args_tuple);
							tmp_term1.second += m_ck2b[j];
							// Create the second term, using the first one's coefficient and the appropriate code.
							Cterm &tmp_term2 = *m_tmp_term2;
							tmp_term2.first = tmp_term1.first;
							tmp_term2.second = m_ck1[i] + m_ck2a[j];
							piranha_assert(tmp_term1.second >= 0);
							piranha_assert(tmp_term2.second >= 0);
							// Now fix flavours and coefficient signs.
							if (f1[i] == f2[j]) {
								if (!f1[i]) {
									tmp_term2.first.invert_sign(m_args_tuple);
								}
								// Insert into cosine container.
								std::pair<bool,c_iterator> res = cms_cos.find(tmp_term1.second);
								if (res.first) {
									res.second->first.add(tmp_term1.first, m_args_tuple);
								} else {
									cms_cos.insert_new(tmp_term1,res.second);
								}
								res = cms_cos.find(tmp_term2.second);
								if (res.first) {
									res.second->first.add(tmp_term2.first, m_args_tuple);
								} else {
									cms_cos.insert_new(tmp_term2,res.second);
								}
							} else {
								if (f1[i]) {
									tmp_term1.first.invert_sign(m_args_tuple);
								}
								// Insert into sine container.
								std::pair<bool,c_iterator> res = cms_sin.find(tmp_term1.second);
								if (res.first) {
									res.second->first.add(tmp_term1.first, m_args_tuple);
								} else {
									cms_sin.insert_new(tmp_term1,res.second);
								}
								res = cms_sin.find(tmp_term2.second);
								if (res.first) {
									res.second->first.add(tmp_term2.first, m_args_tuple);
								} else {
									cms_sin.insert_new(tmp_term2,res.second);
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
						const Ckey			*m_ck2a;
						const Ckey			*m_ck2b;
						const GenericTruncator		&m_trunc;
						std::pair<HashSet *,HashSet *>	*m_cms;
						Cterm				*m_tmp_term1;
						Cterm				*m_tmp_term2;
						const ArgsTuple			&m_args_tuple;
					};
					template <class GenericTruncator>
					void perform_hash_coded_multiplication(const cf_type1 *tc1, const cf_type2 *tc2,
						const term_type1 **t1, const term_type2 **t2, const GenericTruncator &trunc)
					{
						stats::trace_stat("mult_st",std::size_t(0),boost::lambda::_1 + 1);
						typedef coded_hash_table<cf_type1, max_fast_int, std_counting_allocator<char> > csht;
						typedef typename csht::iterator c_iterator;
						// Let's find a sensible size hint.
						const std::size_t n_codes = boost::numeric_cast<std::size_t>(boost::numeric::width(this->m_fast_h) + 1);
						const std::size_t size_hint = static_cast<std::size_t>(
							std::max<double>(this->m_density1,this->m_density2) * n_codes);
						csht cms_cos(size_hint), cms_sin(size_hint);
						std::pair<csht *, csht *> res(&cms_cos,&cms_sin);
						const std::size_t size1 = this->m_terms1.size(), size2 = this->m_terms2.size();
						const args_tuple_type &args_tuple = this->m_args_tuple;
						const max_fast_int *ck1 = &this->m_ckeys1[0], *ck2a = &this->m_ckeys2a[0], *ck2b = &this->m_ckeys2b[0];
						// Find out a suitable block size.
						const std::size_t block_size = this->template compute_block_size<sizeof(std::pair<cf_type1,max_fast_int>)>();
						__PDEBUG(std::cout << "Block size: " << block_size << '\n';)
						std::pair<cf_type1,max_fast_int> tmp_term1, tmp_term2;
						hash_functor<std::pair<cf_type1,max_fast_int>,max_fast_int,GenericTruncator,csht>
							hm(&m_flavours1[0],&m_flavours2[0],tc1,tc2,t1,t2,ck1,ck2a,ck2b,trunc,&res,&tmp_term1,&tmp_term2,args_tuple);
						this->blocked_multiplication(block_size,size1,size2,hm);
						__PDEBUG(std::cout << "Done Poisson series hash coded multiplying\n");
						term_type1 tmp_term;
						{
							const c_iterator c_it_f = cms_cos.end();
							for (c_iterator c_it = cms_cos.begin(); c_it != c_it_f; ++c_it) {
								(c_it->first).divide_by(2,args_tuple);
								this->decode(c_it->first,c_it->second + 2 * this->m_fast_h.lower(),tmp_term);
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
								c_it->first.divide_by(2,args_tuple);
								this->decode(c_it->first,c_it->second + 2 * this->m_fast_h.lower(),tmp_term);
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
