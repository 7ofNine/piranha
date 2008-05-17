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

#ifndef PIRANHA_POLYNOMIAL_H
#define PIRANHA_POLYNOMIAL_H

#include <boost/operators.hpp>
#include <boost/multi_index_container.hpp>
#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include "../base_classes/base_series.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../integer_typedefs.h"
#include "../polynomial_common/common_polynomial_toolbox.h"
#include "../polynomial_common/monomial.h"
#include "../settings.h"

#define POLYNOMIAL_TERM REAL_NAMED_SERIES_TERM(piranha::monomial)
#define POLYNOMIAL REAL_NAMED_SERIES(piranha::polynomial)
#define POLYNOMIAL_BASE_ANCESTOR REAL_NAMED_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial)
#define POLYNOMIAL_NAMED_ANCESTOR REAL_NAMED_SERIES_NAMED_ANCESTOR(boost::tuple<poly_args_descr>,piranha::polynomial)
#define POLYNOMIAL_MULT_ANCESTOR piranha::series_multiplication< POLYNOMIAL, Multiplier, Truncator>
#define POLYNOMIAL_POWER_SERIES_ANCESTOR power_series<0,POLYNOMIAL >
#define POLYNOMIAL_COMMON_POLYNOMIAL_ANCESTOR common_polynomial_toolbox< POLYNOMIAL >
#define POLYNOMIAL_SPECIAL_FUNCTIONS_ANCESTOR named_series_special_functions< POLYNOMIAL >

namespace piranha
{
// NOTE: this is an example of replacing the multiindex container for polynomials with a thin wrapper
// around tr1::unordered_set.
//
//   #include <tr1/unordered_set>
//
//   template <class Term, class Allocator>
//   class polynomial_container_type
//   {
//       typedef std::tr1::unordered_set<Term,typename Term::hasher,std::equal_to<Term>,Allocator> container_type;
//     public:
//       typedef typename container_type::const_iterator const_pinpoint_iterator;
//       typedef typename container_type::iterator pinpoint_iterator;
//       typedef const_pinpoint_iterator const_sorted_iterator;
//       typedef pinpoint_iterator sorted_iterator;
//       polynomial_container_type() {}
//       const_sorted_iterator begin() const {return m_container.begin();}
//       sorted_iterator begin() {return m_container.begin();}
//       const_sorted_iterator end() const {return m_container.end();}
//       sorted_iterator end() {return m_container.end();}
//       void swap(polynomial_container_type &p) {m_container.swap(p.m_container);}
//       size_t size() const {return m_container.size();}
//       bool empty() const {return m_container.empty();}
//       // These cannot be const since we may need to modify the resulting iterator.
//       pinpoint_iterator find(const Term &t) {return m_container.find(t);}
//       sorted_iterator insert(const const_sorted_iterator &, const Term &t)
//       {
//         std::pair<pinpoint_iterator,bool> res(m_container.insert(t));
//         p_assert(res.second);
//         return m_container.begin();
//       }
//       template <class Modifier>
//         bool modify(pinpoint_iterator &it, Modifier &m)
//       {
//         m(*it);
//         return true;
//       }
//       void erase(const const_pinpoint_iterator &it) {m_container.erase(*it);}
//       void max_load_factor(const float &l) {m_container.max_load_factor(l);}
//     private:
//       container_type  m_container;
//   };

	// TODO: generalise here (and elsewhere) the backbone container by introducing thin wrappers around standard containers
	// so that below the typedefs and aliases are truly generic.
	template < NAMED_SERIES_TP_DECL = std::allocator<char> >
	class polynomial:
				public POLYNOMIAL_BASE_ANCESTOR,
				public POLYNOMIAL_NAMED_ANCESTOR,
				public POLYNOMIAL_POWER_SERIES_ANCESTOR,
				public POLYNOMIAL_MULT_ANCESTOR,
				public POLYNOMIAL_COMMON_POLYNOMIAL_ANCESTOR,
				public POLYNOMIAL_SPECIAL_FUNCTIONS_ANCESTOR,
				boost::ring_operators < POLYNOMIAL,
				boost::ring_operators < POLYNOMIAL, max_fast_int,
				boost::ring_operators < POLYNOMIAL, double,
				boost::dividable < POLYNOMIAL, max_fast_int,
				boost::dividable < POLYNOMIAL, double
				> > > > >
	{
			typedef POLYNOMIAL_TERM term_type_;
			typedef Allocator allocator_type;
			typedef POLYNOMIAL_NAMED_ANCESTOR named_ancestor;
			typedef POLYNOMIAL_BASE_ANCESTOR base_ancestor;
			typedef boost::multi_index_container<term_type_, typename I<term_type_>::type, allocator_type> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class POLYNOMIAL_NAMED_ANCESTOR;
			friend class POLYNOMIAL_BASE_ANCESTOR;
			friend class POLYNOMIAL_MULT_ANCESTOR;
			friend class POLYNOMIAL_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class named_series_complex_toolbox<POLYNOMIAL>;
			// Override base_series::real_pow with the one from the common polynomial toolbox.
			using POLYNOMIAL_COMMON_POLYNOMIAL_ANCESTOR::real_pow;
		public:
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef Multiplier<polynomial, polynomial, typename named_ancestor::args_tuple_type, Truncator> multiplier_type;
			// Ctors.
			NAMED_SERIES_CTORS(polynomial);
			// Ctor from psym.
			explicit polynomial(const psym &p) {
				nth_index<1>().max_load_factor(settings::load_factor());
				named_ancestor::template construct_from_psym<0>(p);
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

#define COMPLEX_POLYNOMIAL_TERM COMPLEX_NAMED_SERIES_TERM(piranha::monomial)
#define COMPLEX_POLYNOMIAL COMPLEX_NAMED_SERIES(piranha::polynomial)
#define COMPLEX_POLYNOMIAL_BASE_ANCESTOR COMPLEX_NAMED_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial)
#define COMPLEX_POLYNOMIAL_NAMED_ANCESTOR COMPLEX_NAMED_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::poly_args_descr>, \
		piranha::polynomial)
#define COMPLEX_POLYNOMIAL_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_POLYNOMIAL, Multiplier, Truncator>
#define COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex_toolbox<POLYNOMIAL>

namespace std
{
	template < NAMED_SERIES_TP_DECL >
	class complex<POLYNOMIAL>:
				public COMPLEX_POLYNOMIAL_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_ANCESTOR,
				public COMPLEX_POLYNOMIAL_MULT_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX,
				boost::ring_operators < COMPLEX_POLYNOMIAL,
				boost::ring_operators < COMPLEX_POLYNOMIAL, piranha::max_fast_int,
				boost::ring_operators < COMPLEX_POLYNOMIAL, double,
				boost::dividable < COMPLEX_POLYNOMIAL, piranha::max_fast_int,
				boost::dividable < COMPLEX_POLYNOMIAL, double,
				boost::ring_operators < COMPLEX_POLYNOMIAL, POLYNOMIAL,
				boost::ring_operators < COMPLEX_POLYNOMIAL, complex<piranha::max_fast_int>,
				boost::ring_operators < COMPLEX_POLYNOMIAL, complex<double>,
				boost::dividable < COMPLEX_POLYNOMIAL, complex<piranha::max_fast_int>,
				boost::dividable < COMPLEX_POLYNOMIAL, complex<double>
				> > > > > > > > > >
	{
			typedef COMPLEX_POLYNOMIAL_TERM term_type_;
			typedef Allocator allocator_type;
			typedef COMPLEX_POLYNOMIAL_NAMED_ANCESTOR named_ancestor;
			typedef COMPLEX_POLYNOMIAL_BASE_ANCESTOR base_ancestor;
			typedef boost::multi_index_container<term_type_, typename I<term_type_>::type, allocator_type> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class COMPLEX_POLYNOMIAL_NAMED_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_BASE_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_MULT_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX;
			friend class piranha::base_series_complex_toolbox<POLYNOMIAL>;
		public:
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::add;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::add;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::subtract;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::mult_by;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::divide_by;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator+=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator-=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator*=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator/=;
			// Needed typedefs.
			typedef POLYNOMIAL value_type;
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			// Ctors.
			NAMED_SERIES_CTORS(complex);
			COMPLEX_NAMED_SERIES_CTORS(COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX);
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

// Overload standard math functions for polynomials.
namespace std
{
	template < NAMED_SERIES_TP_DECL >
	POLYNOMIAL pow(const POLYNOMIAL &x, const double &y)
	{
		POLYNOMIAL retval(x.pow(y));
		return retval;
	}
}

#endif
