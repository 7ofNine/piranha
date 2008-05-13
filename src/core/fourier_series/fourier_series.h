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
#include <memory> // For default allocator.

#include "../arg_manager.h"
#include "../base_classes/base_series.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/named_series.h"
#include "../integer_typedefs.h"
#include "../poisson_series_common/poisson_series_term.h"
#include "../ntuple.h"



#define FOURIER_SERIES_TP_DECL class Cf, class Trig, template <class> class I, \
  template <class, class, class, template <class> class> class Multiplier, \
  template <class> class Truncator, class Allocator
#define FOURIER_SERIES_TP Cf,Trig,I,Multiplier,Truncator,Allocator
#define FOURIER_SERIES fourier_series<FOURIER_SERIES_TP>
#define FOURIER_SERIES_BASE_ANCESTOR base_series<poisson_series_term<Cf,Trig,'|',Allocator>,'\n', \
  Allocator,FOURIER_SERIES >
#define FOURIER_SERIES_NAMED_ANCESTOR named_series<boost::tuple<trig_args_descr>,FOURIER_SERIES >
#define FOURIER_SERIES_MULT_ANCESTOR series_multiplication< FOURIER_SERIES, Multiplier, Truncator>

namespace piranha
{
	template < FOURIER_SERIES_TP_DECL = std::allocator<char> >
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
			typedef poisson_series_term < Cf, Trig, '|', Allocator > term_type_;
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
		public:
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			// Ctors.
			__PIRANHA_NAMED_SERIES_CTORS(fourier_series);
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

// __PIRANHA_COMPLEX_FOURIER_SERIES_BASE_ANCESTOR
// 
// namespace std
// {
// 	template < FOURIER_SERIES_TP_DECL >
// 	class complex<piranha::FOURIER_SERIES>:
// 				public FOURIER_SERIES_BASE_ANCESTOR,
// 				public FOURIER_SERIES_NAMED_ANCESTOR,
// 				public FOURIER_SERIES_MULT_ANCESTOR,
// 				boost::ring_operators < FOURIER_SERIES,
// 				boost::ring_operators < FOURIER_SERIES, max_fast_int,
// 				boost::ring_operators < FOURIER_SERIES, double,
// 				boost::dividable < FOURIER_SERIES, max_fast_int,
// 				boost::dividable < FOURIER_SERIES, double
// 				> > > > >
// 		
// }

// Overload standard math functions for Fourier series.
namespace std
{
	// Overload power function for Fourier series.
	template < FOURIER_SERIES_TP_DECL >
	piranha::FOURIER_SERIES pow(const piranha::FOURIER_SERIES &x, const double &y)
	{
		piranha::FOURIER_SERIES retval(x.pow(y));
		return retval;
	}
}

#undef FOURIER_SERIES_TP_DECL
#undef FOURIER_SERIES_TP
#undef FOURIER_SERIES
#undef FOURIER_SERIES_BASE_ANCESTOR
#undef FOURIER_SERIES_NAMED_ANCESTOR
#undef FOURIER_SERIES_MULT_ANCESTOR

#endif
