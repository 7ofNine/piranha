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

#ifndef PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H

#include <algorithm> // For sorting.
#include <boost/type_traits/is_same.hpp>
#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_comparisons.h"
#include "../common_functors.h"
#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../ntuple.h"
#include "../p_assert.h"
#include "../psym.h"
#include "jacobi_anger_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class ArgsTuple>
	class ps_binomial_sorter
	{
		public:
			ps_binomial_sorter(const ArgsTuple &) {}
			template <class Term>
			bool operator()(const Term *t1, const Term *t2) const {
				const max_fast_int d1 = t1->m_cf.min_degree(), d2 = t2->m_cf.min_degree();
				if (d1 == d2) {
					// NOTICE: the idea is that for leading terms with equal
					// min_degree we choose the ones that have
					// unity trig vector, so that we increase the chance of
					// being able to perform the expansion.
					if (t1->m_key.is_unity()) {
						return true;
					} else if (t2->m_key.is_unity()) {
						return false;
					}
					return (t1->m_key < t2->m_key);
				} else {
					return (d1 < d2);
				}
			}
	};

	template <class Derived>
	class common_poisson_series_toolbox:
		public jacobi_anger_toolbox<1, Derived>,
		public binomial_exponentiation_toolbox<Derived,ps_binomial_sorter>
	{
			typedef jacobi_anger_toolbox<1, Derived> jacang_ancestor;
		public:
			// NOTICE: move this into private?
			// NOTICE: this method assumes that the input args tuple already hase merged in as
			// trig arguments the poly arguments (see also below).
			template <class ArgsTuple>
			std::complex<Derived> ei(const ArgsTuple &args_tuple) const {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				typedef typename Derived::term_type term_type;
				typedef typename term_type::cf_type::term_type::cf_type poly_cf_type;
				typedef typename std::complex<Derived>::term_type::cf_type complex_cf_type;
				typedef typename std::vector<term_type const *>::const_iterator const_iterator;
				// Cache and sort the terms.
				std::vector<term_type const *> cache(derived_const_cast->cache_pointers());
				std::sort(cache.begin(),cache.end(),cf_norm_comparison_reverse<ArgsTuple>(args_tuple));
				// Get the term that has unity trig vector and whose coefficient is a linear polynomial with integer
				// coefficients or a linear polynomial with integer coefficients and a single coefficient.
				std::pair<const_iterator, std::pair<std::vector<poly_cf_type>, std::vector<max_fast_int> > >
				int_linear_term(get_int_linear_term<term_type, poly_cf_type>(cache,args_tuple));
				// Expand using Jacobi-Anger's identity.
				std::complex<Derived> retval;
				jacang_ancestor::jacobi_anger(cache, int_linear_term.first, retval, args_tuple);
				// If the linear term was found, take care of it.
				if (int_linear_term.first != cache.end()) {
					std::complex<Derived> tmp_series;
					// Let's build the term to be inserted in tmp_series.
					complex_term_type tmp_term1;
					complex_term_type tmp_term2;
					tmp_term1.m_cf = complex_cf_type(std::complex<max_fast_int>(1, 0), args_tuple);
					tmp_term1.m_key.assign_int_vector(int_linear_term.second.second);
					tmp_term1.m_key.flavour() = true;
					tmp_term2.m_cf = complex_cf_type(std::complex<max_fast_int>(0, 1), args_tuple);
					tmp_term2.m_key.assign_int_vector(int_linear_term.second.second);
					tmp_term2.m_key.flavour() = false;
					tmp_series.insert(tmp_term1, args_tuple);
					tmp_series.insert(tmp_term2, args_tuple);
					// Take care of the numerical-coefficient-only term, if any.
					if (int_linear_term.second.first.size() > 0) {
						std::complex<Derived> tmp_series2;
						// NOTE: tmp_series2's arguments tuple is empty, it does not matter
						// because we are explicitly dealing with single numerical coefficient series.
						complex_term_type tmp_term;
						tmp_term.m_cf.insert(
							typename std::complex<Derived>::term_type::cf_type::term_type(
								int_linear_term.second.first[0].ei(tmp_series2.arguments()),
								typename std::complex<Derived>::term_type::cf_type::term_type::key_type()
							),
							tmp_series2.arguments());
						tmp_series2.insert(tmp_term, tmp_series2.arguments());
						tmp_series.mult_by(tmp_series2, args_tuple);
					}
					retval.mult_by(tmp_series, args_tuple);
				}
				return retval;
			}
			std::complex<Derived> ei() const {
				// In order to account for a potential integer linear combination of arguments
				// we must merge in as trigonometric arguments the polynomial arguments. The safe
				// way to do this is by using named_series::merge_args with a phony series having zero
				// polynomial arguments and as trigonometric arguments the polynomial arguments of this.
				Derived copy(*derived_const_cast), tmp;
				tmp.m_arguments.template get<1>() = derived_const_cast->m_arguments.template get<0>();
				copy.merge_args(tmp);
				// Now we can call in the ei method from above.
				std::complex<Derived> retval(copy.ei(copy.m_arguments));
				retval.m_arguments = copy.arguments();
				retval.trim();
				return retval;
			}
			Derived cos() const {
				return ei().real();
			}
			Derived sin() const {
				return ei().imag();
			}
			// We have to specialise this in order to prepare the arguments tuple for the fact
			// that poly args of s may be added as trig args (as s will be used as argument for sines
			// and/or cosines).
			template <class SubSeries>
			Derived sub(const psym &arg, const SubSeries &series) const {
				typedef typename Derived::args_tuple_type args_tuple_type;
				typedef typename ntuple<std::pair<bool, size_t>, Derived::n_arguments_sets>::type pos_tuple_type;
				typedef typename Derived::term_type::cf_type::
					template sub_cache_selector<SubSeries,typename Derived::term_type::key_type::
					template sub_cache_selector<SubSeries,boost::tuples::null_type,args_tuple_type>
					::type,args_tuple_type>::type sub_caches_type;
				p_static_check(boost::tuples::length<sub_caches_type>::value ==
					boost::tuples::length<pos_tuple_type>::value,
					"Size mismatch for position and cache tuples in Poisson series substitution.");
				Derived this_copy(*derived_const_cast);
				SubSeries s_copy(series), tmp;
				// Assign as tmp's trig arguments series's polynomial arguments.
				args_tuple_type tmp_args;
				tmp_args.template get<1>() = series.arguments().template get<0>();
				tmp.set_arguments(tmp_args);
				// After the next line, s_copy's args layout is compatible with tmp's.
				s_copy.merge_args(tmp);
				// After the next line, this_copy's args layout is compatible with s_copy's
				this_copy.merge_args(s_copy);
				// Finally, have s_copy have compatible arguments with this_copy. This is needed because
				// we will be using this_copy's arguments as args_tuple in all base functions used from now
				// on, including functions taking s_copy as arguments and which do not know anything about
				// this_copy.
				s_copy.merge_args(this_copy);
				// Init sub caches using s and this_copy.m_arguments.
				sub_caches_type sub_caches;
				init_sub_caches<sub_caches_type,SubSeries,args_tuple_type>::run(sub_caches,s_copy,
					&this_copy.m_arguments);
				pos_tuple_type pos_tuple;
				psym_p p(psyms::get_pointer(arg));
				named_series_get_psym_p_positions<pos_tuple_type, args_tuple_type>::run(p, pos_tuple,
					this_copy.m_arguments);
				Derived retval(this_copy.template base_sub<Derived,typename Derived::sub_functor>(pos_tuple, s_copy,
					sub_caches, this_copy.m_arguments));
				retval.m_arguments = this_copy.m_arguments;
				retval.trim();
				return retval;
			}
			template <class SubSeries>
			Derived ei_sub(const psym &arg, const SubSeries &series) const {
				typedef typename Derived::args_tuple_type args_tuple_type;
				typedef typename ntuple<std::pair<bool, size_t>, Derived::n_arguments_sets>::type pos_tuple_type;
				typedef typename Derived::term_type::cf_type::
					template ei_sub_cache_selector<SubSeries,typename Derived::term_type::key_type::
					template ei_sub_cache_selector<SubSeries,boost::tuples::null_type,args_tuple_type>
					::type,args_tuple_type>::type sub_caches_type;
				p_static_check(boost::tuples::length<sub_caches_type>::value ==
					boost::tuples::length<pos_tuple_type>::value,
					"Size mismatch for position and cache tuples in Poisson series ei substitution.");
				Derived this_copy(*derived_const_cast);
				SubSeries s_copy(series);
				this_copy.merge_args(s_copy);
				s_copy.merge_args(this_copy);
				sub_caches_type sub_caches;
				init_sub_caches<sub_caches_type,SubSeries,args_tuple_type>::run(sub_caches,s_copy,
					&this_copy.m_arguments);
				pos_tuple_type pos_tuple;
				psym_p p(psyms::get_pointer(arg));
				named_series_get_psym_p_positions<pos_tuple_type, args_tuple_type>::run(p, pos_tuple,
					this_copy.m_arguments);
				Derived retval(this_copy.template base_sub<Derived,ei_sub_functor>(pos_tuple, s_copy,
					sub_caches, this_copy.m_arguments));
				retval.m_arguments = this_copy.m_arguments;
				retval.trim();
				return retval;
			}
			template <class FourierSeries>
			FourierSeries to_fs() const {
				typedef typename Derived::const_iterator const_iterator;
				typedef typename FourierSeries::term_type fourier_term;
				typename ntuple<vector_psym_p,1>::type args_tuple(derived_const_cast->arguments().template get<1>());
				FourierSeries retval;
				retval.set_arguments(args_tuple);
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) {
					if (!it->m_cf.is_single_cf()) {
						throw unsuitable("Polynomial coefficient cannot be converted to numerical coefficient.");
					}
					retval.insert(fourier_term(typename fourier_term::cf_type(it->m_cf.begin()->m_cf),
						typename fourier_term::key_type(it->m_key)),args_tuple);
				}
				return retval;
			}
		private:
			template <class Term, class PolyCf, class ArgsTuple>
			std::pair<typename std::vector<Term const *>::const_iterator, std::pair<std::vector<PolyCf>,
			std::vector<max_fast_int> > >
			static get_int_linear_term(const std::vector<const Term *> &v, const ArgsTuple &args_tuple) {
				p_static_check((boost::is_same<PolyCf,
					typename Derived::term_type::cf_type::term_type::cf_type>::value),
					"Coefficient type mismatch in Poisson series toolbox.");
				typedef typename std::vector<Term const *>::const_iterator const_iterator;
				const const_iterator it_f = v.end();
				std::pair<const_iterator, std::pair<std::vector<PolyCf>, std::vector<max_fast_int> > > retval;
				retval.first = it_f;
				// Make space to accomodate all the elements of the linear combination.
				// We need as much space as the number of trig args.
				retval.second.second.resize(args_tuple.template get<1>().size());
				for (const_iterator it = v.begin(); it != it_f; ++it) {
					// If the term's trigonometric part is unity, let's see if we can extract a
					// linear combination of arguments from the corresponding polynomial.
					if ((*it)->m_key.is_unity()) {
						try {
							(*it)->m_cf.template get_int_linear_combination<1>(retval.second, args_tuple);
						} catch (const unsuitable &) {
							// If we are unable to extract a proper linear combination from the unity term, erase retval
							// and break out.
							retval.second.first.clear();
							retval.second.second.clear();
							break;
						}
						retval.first = it;
						break;
					}
				}
				p_assert(retval.second.first.size() <= 1);
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
