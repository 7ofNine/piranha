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
#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/common_indices.h"
#include "../base_classes/series_multiindex_backend.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../integer_typedefs.h"
#include "common_fourier_series_toolbox.h"
#include "fourier_series_term.h"

#define FOURIER_SERIES_TERM E0_SERIES_TERM(piranha::fourier_series_term)
#define FOURIER_SERIES E0_SERIES(piranha::fourier_series)
#define FOURIER_SERIES_BASE_ANCESTOR E0_SERIES_BASE_ANCESTOR(piranha::fourier_series_term,piranha::fourier_series)
#define FOURIER_SERIES_NAMED_ANCESTOR E0_SERIES_NAMED_ANCESTOR(boost::tuple<trig_args_descr>,piranha::fourier_series)
#define FOURIER_SERIES_MULT_ANCESTOR piranha::series_multiplication< FOURIER_SERIES, Multiplier, Truncator>
#define FOURIER_SERIES_COMMON_ANCESTOR piranha::common_fourier_series_toolbox< FOURIER_SERIES >
#define FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR piranha::binomial_exponentiation_toolbox< FOURIER_SERIES >
#define FOURIER_SERIES_SPECIAL_FUNCTIONS_ANCESTOR piranha::named_series_special_functions< FOURIER_SERIES >

namespace piranha
{
	template < E0_SERIES_TP_DECL = std::allocator<char> >
	class fourier_series:
				public FOURIER_SERIES_BASE_ANCESTOR,
				public FOURIER_SERIES_NAMED_ANCESTOR,
				public FOURIER_SERIES_MULT_ANCESTOR,
				public FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR,
				public FOURIER_SERIES_COMMON_ANCESTOR,
				public FOURIER_SERIES_SPECIAL_FUNCTIONS_ANCESTOR,
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
			typedef series_multiindex_backend < term_type_, random_keyhash_index,
			typename allocator_type::template rebind<term_type_>::other > container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class FOURIER_SERIES_NAMED_ANCESTOR;
			friend class FOURIER_SERIES_BASE_ANCESTOR;
			friend class FOURIER_SERIES_MULT_ANCESTOR;
			friend class FOURIER_SERIES_COMMON_ANCESTOR;
			friend class FOURIER_SERIES_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class named_series_complex_toolbox<FOURIER_SERIES>;
			using FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::real_power;
			using FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::negative_integer_power;
			using FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::nth_root;
		public:
			// TODO: do we really need _all_ those double interfaces here? Some can probably be dropped.
			using named_ancestor::norm;
			using base_ancestor::norm;
			using named_ancestor::eval;
			using base_ancestor::eval;
			using named_ancestor::pow;
			using base_ancestor::pow;
			using named_ancestor::root;
			using base_ancestor::root;
			using named_ancestor::partial;
			using base_ancestor::partial;
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef typename Multiplier::template get_type < fourier_series, fourier_series,
			args_tuple_type, Truncator > multiplier_type;
			// Ctors.
			NAMED_SERIES_CTORS(fourier_series, 0);
			// Index interface.
			SERIES_INDEX_INTERFACE;
		private:
			container_type  m_container;
	};
}

#define COMPLEX_FOURIER_SERIES_TERM COMPLEX_E0_SERIES_TERM(piranha::fourier_series_term)
#define COMPLEX_FOURIER_SERIES COMPLEX_E0_SERIES(piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_BASE_ANCESTOR COMPLEX_E0_SERIES_BASE_ANCESTOR(piranha::fourier_series_term,piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR COMPLEX_E0_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::trig_args_descr>, \
		piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_FOURIER_SERIES, Multiplier, Truncator>
#define COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX piranha::base_series_complex_toolbox<FOURIER_SERIES>
#define COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex_toolbox<FOURIER_SERIES>
#define COMPLEX_FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR piranha::binomial_exponentiation_toolbox< COMPLEX_FOURIER_SERIES >
#define COMPLEX_FOURIER_SERIES_SPECIAL_FUNCTIONS_ANCESTOR piranha::named_series_special_functions< COMPLEX_FOURIER_SERIES >

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<FOURIER_SERIES>:
				public COMPLEX_FOURIER_SERIES_BASE_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_MULT_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX,
				public COMPLEX_FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_SPECIAL_FUNCTIONS_ANCESTOR,
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
			typedef piranha::series_multiindex_backend < term_type_, piranha::random_keyhash_index,
			typename allocator_type::template rebind<term_type_>::other > container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_BASE_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_MULT_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX;
			friend class COMPLEX_FOURIER_SERIES_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX;
			friend class piranha::common_fourier_series_toolbox<FOURIER_SERIES>;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::real_power;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::negative_integer_power;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::nth_root;
		public:
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::real;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::real;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::imag;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::imag;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::construct_from_real;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::construct_from_real;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::construct_from_real_imag;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::construct_from_real_imag;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::add;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::add;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::subtract;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::mult_by;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::divide_by;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator+=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator-=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator*=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator/=;
			using named_ancestor::norm;
			using base_ancestor::norm;
			using named_ancestor::eval;
			using base_ancestor::eval;
			using named_ancestor::pow;
			using base_ancestor::pow;
			using named_ancestor::root;
			using base_ancestor::root;
			using named_ancestor::partial;
			using base_ancestor::partial;
			// Needed typedefs.
			typedef FOURIER_SERIES value_type;
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef typename Multiplier::template get_type < complex, complex,
			args_tuple_type, Truncator > multiplier_type;
			NAMED_SERIES_CTORS(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX);
			SERIES_INDEX_INTERFACE;
		private:
			container_type  m_container;
	};
}

// Overload standard math functions for Fourier series.
namespace std
{
	// Overload power function for Fourier series.
	template < E0_SERIES_TP_DECL >
	FOURIER_SERIES pow(const FOURIER_SERIES &x, const double &y)
	{
		FOURIER_SERIES retval(x.pow(y));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	FOURIER_SERIES pow(const FOURIER_SERIES &x, const piranha::max_fast_int &n)
	{
		FOURIER_SERIES retval(x.pow(n));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	COMPLEX_FOURIER_SERIES pow(const COMPLEX_FOURIER_SERIES &x, const double &y)
	{
		COMPLEX_FOURIER_SERIES retval(x.pow(y));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	COMPLEX_FOURIER_SERIES pow(const COMPLEX_FOURIER_SERIES &x, const piranha::max_fast_int &n)
	{
		COMPLEX_FOURIER_SERIES retval(x.pow(n));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	FOURIER_SERIES cos(const FOURIER_SERIES &x)
	{
		return x.cos();
	}

	template < E0_SERIES_TP_DECL >
	FOURIER_SERIES sin(const FOURIER_SERIES &x)
	{
		return x.sin();
	}
}

#endif
