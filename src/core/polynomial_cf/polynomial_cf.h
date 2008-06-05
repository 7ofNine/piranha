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

#ifndef PIRANHA_POLYNOMIAL_CF_H
#define PIRANHA_POLYNOMIAL_CF_H

#include <complex>
#include <utility>

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/cf_series.h"
#include "../base_classes/power_series.h"
#include "../base_classes/series_multiindex_backend.h"
#include "../base_classes/series_multiplication.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../polynomial_common/monomial.h"
#include "../proxies.h"

#define POLYNOMIAL_CF_TERM CF_SERIES_TERM(piranha::monomial,'!')
#define POLYNOMIAL_CF E0_SERIES(piranha::polynomial_cf)
#define POLYNOMIAL_CF_BASE_ANCESTOR CF_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial_cf,'!',',')
#define POLYNOMIAL_CF_CF_ANCESTOR piranha::cf_series< POLYNOMIAL_CF >
#define POLYNOMIAL_CF_MULT_ANCESTOR piranha::series_multiplication< POLYNOMIAL_CF, Multiplier, Truncator>
#define POLYNOMIAL_CF_POWER_SERIES_ANCESTOR power_series<0,POLYNOMIAL_CF >
#define POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR binomial_exponentiation_toolbox< POLYNOMIAL_CF >
#define POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR piranha::base_series_special_functions< POLYNOMIAL_CF >

namespace piranha
{
	// NOTE: this assumes that exponents are in position 0 of arguments tuple.
	template <E0_SERIES_TP_DECL>
	class polynomial_cf:
				public POLYNOMIAL_CF_BASE_ANCESTOR,
				public POLYNOMIAL_CF_CF_ANCESTOR,
				public POLYNOMIAL_CF_POWER_SERIES_ANCESTOR,
				public POLYNOMIAL_CF_MULT_ANCESTOR,
				public POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR,
				public POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR
	{
			typedef POLYNOMIAL_CF_TERM term_type_;
			typedef typename term_type_::cf_type cf_type;
			typedef typename term_type_::key_type key_type;
			typedef Allocator allocator_type;
			typedef POLYNOMIAL_CF_CF_ANCESTOR cf_ancestor;
			typedef POLYNOMIAL_CF_BASE_ANCESTOR base_ancestor;
			typedef series_multiindex_backend<term_type_, I, allocator_type> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			friend class POLYNOMIAL_CF_CF_ANCESTOR;
			friend class POLYNOMIAL_CF_BASE_ANCESTOR;
			friend class POLYNOMIAL_CF_MULT_ANCESTOR;
			// Specify we will use the real_power from the binomial exponentiation toolbox.
			using POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR::real_power;
			using POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR::negative_integer_power;
			using POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR::nth_root;
		public:
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef Multiplier<polynomial_cf, polynomial_cf, boost::tuples::null_type, Truncator> multiplier_type;
			CF_SERIES_CTORS(polynomial_cf);
			template <class ArgsTuple>
			explicit polynomial_cf(const psym_p &p, const int &n, const ArgsTuple &a) {
				base_ancestor::construct_from_psym_p(p, n, a);
			}
			// Needed getters and setters.
			template <int N>
			typename container_type::template nth_index<N>::type &nth_index() {
				return m_container.template get<N>();
			}
			template <int N>
			const typename container_type::template nth_index<N>::type &nth_index() const {
				return m_container.template get<N>();
			}
			// TODO: place some of these methods into common polynomial toolbox?
			/// Return a single coefficient and a vector of integers representing the polynomial.
			void get_int_linear_combination(std::pair<std::vector<cf_type>, std::vector<max_fast_int> > &res) const {
				const const_sorted_iterator it_f = nth_index<0>().end();
				for (const_sorted_iterator it = nth_index<0>().begin(); it != it_f; ++it) {
					const max_fast_int pos = it->m_key.linear_arg_position();
					if (pos >= 0) {
						p_assert(pos < (max_fast_int)res.second.size());
						res.second[(size_t)pos] = it->m_cf.get_int();
					} else {
						res.first.push_back(it->m_cf);
					}
				}
			}
		private:
			container_type  m_container;
	};
}

#define COMPLEX_POLYNOMIAL_CF_TERM COMPLEX_CF_SERIES_TERM(piranha::monomial,'!')
#define COMPLEX_POLYNOMIAL_CF COMPLEX_E0_SERIES(piranha::polynomial_cf)
#define COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR COMPLEX_CF_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial_cf,'!',',')
#define COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR piranha::cf_series< COMPLEX_POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_POLYNOMIAL_CF, Multiplier, Truncator>
#define COMPLEX_POLYNOMIAL_CF_POWER_SERIES_ANCESTOR piranha::power_series<0,COMPLEX_POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX piranha::base_series_complex_toolbox< POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR piranha::base_series_special_functions< COMPLEX_POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR piranha::binomial_exponentiation_toolbox< COMPLEX_POLYNOMIAL_CF >

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL_CF>:
				public COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_POWER_SERIES_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_MULT_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR
	{
			typedef COMPLEX_POLYNOMIAL_CF_TERM term_type_;
			typedef typename term_type_::cf_type cf_type;
			typedef typename term_type_::key_type key_type;
			typedef Allocator allocator_type;
			typedef COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR cf_ancestor;
			typedef COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR base_ancestor;
			typedef piranha::series_multiindex_backend<term_type_, I, allocator_type> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			friend class COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_CF_MULT_ANCESTOR;
			friend class piranha::base_series_complex_toolbox<POLYNOMIAL_CF>;
			// Specify we will use the real_power from the binomial exponentiation toolbox.
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR::real_power;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR::negative_integer_power;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_EXPONENTIATION_ANCESTOR::nth_root;
		public:
			using COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX::add;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::add;
			using COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::subtract;
			using COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::mult_by;
			using COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::divide_by;
			// Needed typedefs.
			typedef POLYNOMIAL_CF value_type;
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef Multiplier<complex, complex, boost::tuples::null_type, Truncator> multiplier_type;
			CF_SERIES_CTORS(complex);
			COMPLEX_CF_SERIES_CTORS(COMPLEX_POLYNOMIAL_CF_COMPLEX_TOOLBOX);
			template <class ArgsTuple>
			explicit complex(const piranha::psym_p &p, const int &n, const ArgsTuple &a) {
				base_ancestor::construct_from_psym_p(p, n, a);
			}
			// Needed getters and setters.
			template <int N>
			typename container_type::template nth_index<N>::type &nth_index() {
				return m_container.template get<N>();
			}
			template <int N>
			const typename container_type::template nth_index<N>::type &nth_index() const {
				return m_container.template get<N>();
			}
		private:
			container_type  m_container;
	};
}

namespace piranha
{
	// Specialisation of cf mult proxy for polynomial_cf to use reference.
	template < E0_SERIES_TP_DECL >
	class cf_mult_proxy< POLYNOMIAL_CF >:
				public reference_proxy< POLYNOMIAL_CF >
	{
			typedef reference_proxy< POLYNOMIAL_CF > ancestor;
		public:
			cf_mult_proxy(): ancestor() {}
			void operator=(const POLYNOMIAL_CF &cf) {
				ancestor::assignment(cf);
			}
	};

	// Specialisation of cf mult proxy for complex polynomial_cf to use reference.
	template < E0_SERIES_TP_DECL >
	class cf_mult_proxy< COMPLEX_POLYNOMIAL_CF >:
				public reference_proxy< COMPLEX_POLYNOMIAL_CF >
	{
			typedef reference_proxy< COMPLEX_POLYNOMIAL_CF > ancestor;
		public:
			cf_mult_proxy(): ancestor() {}
			void operator=(const COMPLEX_POLYNOMIAL_CF &cf) {
				ancestor::assignment(cf);
			}
	};
}


#endif
