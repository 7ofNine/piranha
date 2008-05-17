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

#ifndef PIRANHA_FOURIER_SERIES_H
#define PIRANHA_FOURIER_SERIES_H

#include <boost/operators.hpp>
#include <boost/multi_index_container.hpp>
#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include "../arg_manager.h"
#include "../base_classes/base_series.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../integer_typedefs.h"
#include "../poisson_series_common/poisson_series_term.h"
#include "../ntuple.h"

#define FOURIER_SERIES_TERM REAL_NAMED_SERIES_TERM(piranha::poisson_series_term)
#define FOURIER_SERIES REAL_NAMED_SERIES(piranha::fourier_series)
#define FOURIER_SERIES_BASE_ANCESTOR REAL_NAMED_SERIES_BASE_ANCESTOR(piranha::poisson_series_term,piranha::fourier_series)
#define FOURIER_SERIES_NAMED_ANCESTOR REAL_NAMED_SERIES_NAMED_ANCESTOR(boost::tuple<trig_args_descr>,piranha::fourier_series)
#define FOURIER_SERIES_MULT_ANCESTOR piranha::series_multiplication< FOURIER_SERIES, Multiplier, Truncator>

namespace piranha
{
	template < NAMED_SERIES_TP_DECL = std::allocator<char> >
	class fourier_series:
				public FOURIER_SERIES_BASE_ANCESTOR,
				public FOURIER_SERIES_NAMED_ANCESTOR,
				public FOURIER_SERIES_MULT_ANCESTOR,
				boost::ring_operators < FOURIER_SERIES,
				boost::ring_operators < FOURIER_SERIES, max_fast_int,
				boost::ring_operators < FOURIER_SERIES, double,
				boost::dividable < FOURIER_SERIES, max_fast_int,
				boost::dividable < FOURIER_SERIES, double
				> > > > >
	{
			typedef FOURIER_SERIES_TERM term_type_;
			typedef Allocator allocator_type;
			typedef FOURIER_SERIES_NAMED_ANCESTOR named_ancestor;
			typedef FOURIER_SERIES_BASE_ANCESTOR base_ancestor;
			typedef boost::multi_index_container<term_type_, typename I<term_type_>::type, allocator_type> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class FOURIER_SERIES_NAMED_ANCESTOR;
			friend class FOURIER_SERIES_BASE_ANCESTOR;
			friend class FOURIER_SERIES_MULT_ANCESTOR;
			friend class named_series_complex_toolbox<FOURIER_SERIES>;
		public:
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			// Ctors.
			NAMED_SERIES_CTORS(fourier_series);
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

#define COMPLEX_FOURIER_SERIES_TERM COMPLEX_NAMED_SERIES_TERM(piranha::poisson_series_term)
#define COMPLEX_FOURIER_SERIES COMPLEX_NAMED_SERIES(piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_BASE_ANCESTOR COMPLEX_NAMED_SERIES_BASE_ANCESTOR(piranha::poisson_series_term,piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR COMPLEX_NAMED_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::trig_args_descr>, \
		piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_FOURIER_SERIES, Multiplier, Truncator>
#define COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex_toolbox<FOURIER_SERIES>

namespace std
{
	template < NAMED_SERIES_TP_DECL >
	class complex<FOURIER_SERIES>:
				public COMPLEX_FOURIER_SERIES_BASE_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_MULT_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX,
				boost::ring_operators < COMPLEX_FOURIER_SERIES,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, piranha::max_fast_int,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, double,
				boost::dividable < COMPLEX_FOURIER_SERIES, piranha::max_fast_int,
				boost::dividable < COMPLEX_FOURIER_SERIES, double,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, FOURIER_SERIES,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, complex<piranha::max_fast_int>,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, complex<double>,
				boost::dividable < COMPLEX_FOURIER_SERIES, complex<piranha::max_fast_int>,
				boost::dividable < COMPLEX_FOURIER_SERIES, complex<double>
				> > > > > > > > > >
	{
			typedef COMPLEX_FOURIER_SERIES_TERM term_type_;
			typedef Allocator allocator_type;
			typedef COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR named_ancestor;
			typedef COMPLEX_FOURIER_SERIES_BASE_ANCESTOR base_ancestor;
			typedef boost::multi_index_container<term_type_, typename I<term_type_>::type, allocator_type> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_BASE_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_MULT_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX;
			friend class piranha::base_series_complex_toolbox<FOURIER_SERIES>;
		public:
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::add;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::add;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::subtract;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::mult_by;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::divide_by;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator+=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator-=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator*=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator/=;
			// Needed typedefs.
			typedef FOURIER_SERIES value_type;
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			// Ctors.
			NAMED_SERIES_CTORS(complex);
			COMPLEX_NAMED_SERIES_CTORS(COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX);
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

// Overload standard math functions for Fourier series.
namespace std
{
	// Overload power function for Fourier series.
	template < NAMED_SERIES_TP_DECL >
	FOURIER_SERIES pow(const FOURIER_SERIES &x, const double &y)
	{
		FOURIER_SERIES retval(x.pow(y));
		return retval;
	}
}

#endif
