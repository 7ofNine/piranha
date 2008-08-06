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

#ifndef PIRANHA_COMMON_FOURIER_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_FOURIER_SERIES_TOOLBOX_H

#include <algorithm> // For sorting.
#include <complex>
#include <vector>

#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_comparisons.h"
#include "../poisson_series_common/jacobi_anger_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class ArgsTuple>
	class fs_binomial_sorter
	{
		public:
			fs_binomial_sorter(const ArgsTuple &args_tuple):m_args_tuple(args_tuple) {}
			template <class Term>
			bool operator()(const Term &t1, const Term &t2) const {
				const double n1(t1.m_cf.norm(m_args_tuple)), n2(t2.m_cf.norm(m_args_tuple));
				if (n1 != n2) {
					return (n1 > n2);
				} else {
					if (t1.m_key.is_unity()) {
						return true;
					} else if (t2.m_key.is_unity()) {
						return false;
					}
					return (t1.m_key < t2.m_key);
				}
			}
		private:
			const ArgsTuple &m_args_tuple;
	};

	template <class Derived>
	class common_fourier_series_toolbox:
		public jacobi_anger_toolbox<0, Derived>,
		public binomial_exponentiation_toolbox<Derived,fs_binomial_sorter>
	{
			typedef jacobi_anger_toolbox<0, Derived> jacang_ancestor;
		public:
			std::complex<Derived> complexp() const {
				std::complex<Derived> retval(complexp(derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			template <class ArgsTuple>
			std::complex<Derived> complexp(const ArgsTuple &args_tuple) const {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				typedef typename complex_term_type::key_type key_type;
				typedef typename Derived::term_type term_type;
				typedef typename Derived::term_proxy_type term_proxy_type;
				typedef typename std::vector<term_proxy_type>::const_iterator const_iterator;
				std::complex<Derived> retval;
				if (derived_const_cast->is_single_cf()) {
					retval.insert(complex_term_type(derived_const_cast->begin()->
													m_cf.complexp(args_tuple), key_type()),
								  args_tuple);
				} else {
					// Cache and sort the term proxies list. Sorting is reverse (small --> big norms) because
					// in jacang is better to do the small terms first.
					std::vector<term_proxy_type> cache(derived_const_cast->cache_proxies());
					std::sort(cache.begin(),cache.end(),cf_norm_comparison_reverse<ArgsTuple>(args_tuple));
					// Let's find out if there is a constant term. If there is one, it will be skipped
					// and multiplied by the result of the Jacobi-Anger expansion of the other terms later.
					// We treat it this way because the constant term may be a phase with arbitrary value,
					// whereas jacang works for coefficients < 1 or so. Since a constant term has all trig
					// multipliers equal to zero, we can just compute its complex exponential without worrying
					// about convergence.
					const_iterator it = cache.begin();
					const const_iterator it_f = cache.end();
					for (; it != it_f; ++it) {
						if (it->m_key.is_unity()) {
							break;
						}
					}
					// Expand using Jacobi-Anger's identity.
					jacang_ancestor::jacobi_anger(cache, it, retval, args_tuple);
					if (it != it_f) {
						// Take care of the constant element.
						std::complex<Derived> tmp;
						tmp.insert(complex_term_type(it->m_cf.complexp(args_tuple),key_type()),args_tuple);
						retval.mult_by(tmp,args_tuple);
					}
				}
				return retval;
			}
			Derived cos() const {
				return complexp().real();
			}
			Derived sin() const {
				return complexp().imag();
			}
			static std::complex<Derived> Ynm(const max_fast_int &n, const max_fast_int &m,
				const Derived &theta_, const Derived &phi) {
				Derived theta(theta_);
				theta.merge_args(phi);
				std::complex<Derived> retval(Derived::Ynm(n,m,theta,phi,theta.arguments()));
				retval.m_arguments = theta.arguments();
				retval.trim();
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
