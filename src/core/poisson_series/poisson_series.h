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

#ifndef PIRANHA_POISSON_SERIES_H
#define PIRANHA_POISSON_SERIES_H

#include <boost/operators.hpp>
#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/common_indices.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/power_series.h"
#include "../base_classes/series_multiindex_backend.h"
#include "../base_classes/series_multiplication.h"
#include "../integer_typedefs.h"
#include "../fourier_series/fourier_series_term.h"
#include "../poisson_series_common/common_poisson_series_toolbox.h"
#include "../poisson_series_common/celmec_toolbox.h"
#include "../polynomial_cf/polynomial_cf.h"

#define POISSON_SERIES E1_SERIES(piranha::poisson_series)
#define POISSON_SERIES_POLYNOMIAL_CF E1_SERIES_COEFFICIENT(piranha::polynomial_cf)
#define POISSON_SERIES_TERM E1_SERIES_TERM(piranha::fourier_series_term,POISSON_SERIES_POLYNOMIAL_CF)
#define POISSON_SERIES_BASE_ANCESTOR E1_SERIES_BASE_ANCESTOR(piranha::fourier_series_term,POISSON_SERIES_POLYNOMIAL_CF,POISSON_SERIES)
#define POISSON_SERIES_NAMED_ANCESTOR E1_SERIES_NAMED_ANCESTOR(piranha::poly_args_descr,piranha::trig_args_descr,POISSON_SERIES)
#define POISSON_SERIES_MULT_ANCESTOR piranha::series_multiplication< POISSON_SERIES, Mult1, Trunc1>
#define POISSON_SERIES_COMMON_ANCESTOR piranha::common_poisson_series_toolbox< POISSON_SERIES >
#define POISSON_SERIES_POWER_SERIES_ANCESTOR piranha::power_series<0, 0, POISSON_SERIES >
#define POISSON_SERIES_SPECIAL_FUNCTIONS_ANCESTOR piranha::named_series_special_functions< POISSON_SERIES >
#define POISSON_SERIES_CELMEC_ANCESTOR celmec_toolbox< POISSON_SERIES >
#define POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR binomial_exponentiation_toolbox< POISSON_SERIES >

namespace piranha
{
	template < E1_SERIES_TP_DECL = std::allocator<char> >
	class poisson_series:
				public POISSON_SERIES_BASE_ANCESTOR,
				public POISSON_SERIES_NAMED_ANCESTOR,
				public POISSON_SERIES_MULT_ANCESTOR,
				public POISSON_SERIES_COMMON_ANCESTOR,
				public POISSON_SERIES_POWER_SERIES_ANCESTOR,
				public POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR,
				public POISSON_SERIES_SPECIAL_FUNCTIONS_ANCESTOR,
				public POISSON_SERIES_CELMEC_ANCESTOR,
				boost::ring_operators < POISSON_SERIES,
				boost::ring_operators < POISSON_SERIES, max_fast_int,
				boost::ring_operators < POISSON_SERIES, double,
				boost::dividable < POISSON_SERIES, max_fast_int,
				boost::dividable < POISSON_SERIES, double
				> > > > >
	{
			typedef POISSON_SERIES_TERM term_type_;
			typedef typename term_type_::cf_type cf_type;
			typedef typename term_type_::key_type key_type;
			typedef Allocator allocator_type;
			typedef POISSON_SERIES_NAMED_ANCESTOR named_ancestor;
			typedef POISSON_SERIES_BASE_ANCESTOR base_ancestor;
			typedef piranha::series_multiindex_backend<term_type_, piranha::random_keyhash_index,
				typename allocator_type::template rebind<term_type_>::other> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class POISSON_SERIES_NAMED_ANCESTOR;
			friend class POISSON_SERIES_BASE_ANCESTOR;
			friend class POISSON_SERIES_MULT_ANCESTOR;
			friend class POISSON_SERIES_COMMON_ANCESTOR;
			friend class POISSON_SERIES_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class POISSON_SERIES_CELMEC_ANCESTOR;
			friend class named_series_complex_toolbox< POISSON_SERIES >;
			using POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::real_power;
			using POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::negative_integer_power;
			using POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::nth_root;
		public:
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
			using POISSON_SERIES_COMMON_ANCESTOR::sub;
			using base_ancestor::sub;
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef typename Mult1::template get_type<poisson_series, poisson_series, typename named_ancestor::args_tuple_type, Trunc1> multiplier_type;
			// Ctors.
			NAMED_SERIES_CTORS(poisson_series);
			// Ctor from psym.
			explicit poisson_series(const psym &p) {
				named_ancestor::template construct_from_psym<0>(p);
			}
			SERIES_INDEX_INTERFACE;
		private:
			container_type  m_container;
	};
}

#define COMPLEX_POISSON_SERIES std::complex<POISSON_SERIES>
#define COMPLEX_POISSON_SERIES_POLYNOMIAL_CF piranha::polynomial_cf<Cf,Key0,Mult0,Trunc0,Allocator>
#define COMPLEX_POISSON_SERIES_TERM piranha::fourier_series_term<std::complex<COMPLEX_POISSON_SERIES_POLYNOMIAL_CF>,Key1,'|',Allocator>
#define COMPLEX_POISSON_SERIES_BASE_ANCESTOR piranha::base_series<COMPLEX_POISSON_SERIES_TERM,'\n',Allocator,COMPLEX_POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_NAMED_ANCESTOR piranha::named_series<boost::tuple<piranha::poly_args_descr,piranha::trig_args_descr>, \
	COMPLEX_POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_POISSON_SERIES, Mult1, Trunc1>
#define COMPLEX_POISSON_SERIES_COMMON_ANCESTOR piranha::common_poisson_series_toolbox< COMPLEX_POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_POWER_SERIES_ANCESTOR piranha::power_series<0, 0, COMPLEX_POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX piranha::base_series_complex_toolbox< POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex_toolbox< POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR piranha::binomial_exponentiation_toolbox< COMPLEX_POISSON_SERIES >
#define COMPLEX_POISSON_SERIES_SPECIAL_FUNCTIONS_ANCESTOR piranha::named_series_special_functions< POISSON_SERIES >

namespace std
{
	template <E1_SERIES_TP_DECL>
	class complex<POISSON_SERIES>:
				public COMPLEX_POISSON_SERIES_BASE_ANCESTOR,
				public COMPLEX_POISSON_SERIES_NAMED_ANCESTOR,
				public COMPLEX_POISSON_SERIES_MULT_ANCESTOR,
				public COMPLEX_POISSON_SERIES_COMMON_ANCESTOR,
				public COMPLEX_POISSON_SERIES_POWER_SERIES_ANCESTOR,
				public COMPLEX_POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR,
				public COMPLEX_POISSON_SERIES_SPECIAL_FUNCTIONS_ANCESTOR,
				public COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX,
				boost::ring_operators < COMPLEX_POISSON_SERIES,
				boost::ring_operators < COMPLEX_POISSON_SERIES, piranha::max_fast_int,
				boost::ring_operators < COMPLEX_POISSON_SERIES, double,
				boost::dividable < COMPLEX_POISSON_SERIES, piranha::max_fast_int,
				boost::dividable < COMPLEX_POISSON_SERIES, double,
				boost::ring_operators < COMPLEX_POISSON_SERIES, POISSON_SERIES,
				boost::ring_operators < COMPLEX_POISSON_SERIES, complex<piranha::max_fast_int>,
				boost::ring_operators < COMPLEX_POISSON_SERIES, complex<double>,
				boost::dividable < COMPLEX_POISSON_SERIES, complex<piranha::max_fast_int>,
				boost::dividable < COMPLEX_POISSON_SERIES, complex<double>
				> > > > > > > > > >
	{
			typedef COMPLEX_POISSON_SERIES_TERM term_type_;
			typedef typename term_type_::cf_type cf_type;
			typedef typename term_type_::key_type key_type;
			typedef Allocator allocator_type;
			typedef COMPLEX_POISSON_SERIES_NAMED_ANCESTOR named_ancestor;
			typedef COMPLEX_POISSON_SERIES_BASE_ANCESTOR base_ancestor;
			typedef piranha::series_multiindex_backend<term_type_, piranha::random_keyhash_index,
				typename allocator_type::template rebind<term_type_>::other> container_type;
			typedef typename container_type::template nth_index<0>::type sorted_index;
			typedef typename container_type::template nth_index<1>::type pinpoint_index;
			typedef typename named_ancestor::args_tuple_type args_tuple_type;
			friend class COMPLEX_POISSON_SERIES_NAMED_ANCESTOR;
			friend class COMPLEX_POISSON_SERIES_BASE_ANCESTOR;
			friend class COMPLEX_POISSON_SERIES_MULT_ANCESTOR;
			friend class COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX;
			friend class COMPLEX_POISSON_SERIES_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX;
			friend class POISSON_SERIES_COMMON_ANCESTOR;
			using COMPLEX_POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::real_power;
			using COMPLEX_POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::negative_integer_power;
			using COMPLEX_POISSON_SERIES_BINOMIAL_EXPONENTIATION_ANCESTOR::nth_root;
		public:
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::real;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::real;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::imag;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::imag;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::construct_from_real;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::construct_from_real;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::construct_from_real_imag;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::construct_from_real_imag;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::add;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::add;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::subtract;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::mult_by;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::divide_by;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator+=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator-=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator*=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator/=;
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
			using COMPLEX_POISSON_SERIES_COMMON_ANCESTOR::sub;
			using base_ancestor::sub;
			// Needed typedefs.
			typedef POISSON_SERIES value_type;
			// Needed typedefs.
			typedef term_type_ term_type;
			typedef typename sorted_index::const_iterator const_sorted_iterator;
			typedef typename sorted_index::iterator sorted_iterator;
			typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
			typedef typename pinpoint_index::iterator pinpoint_iterator;
			typedef typename Mult1::template get_type<complex, complex, typename named_ancestor::args_tuple_type, Trunc1> multiplier_type;
			// Ctors.
			NAMED_SERIES_CTORS(complex);
			COMPLEX_NAMED_SERIES_CTORS(COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX);
			// Ctor from psym.
			explicit complex(const piranha::psym &p) {
				named_ancestor::template construct_from_psym<0>(p);
			}
			SERIES_INDEX_INTERFACE;
		private:
			container_type  m_container;
	};
}

// Overload standard math functions for Poisson series.
namespace std
{
	template <E1_SERIES_TP_DECL>
	POISSON_SERIES cos(const POISSON_SERIES &p)
	{
		POISSON_SERIES retval = p.cos();
		return retval;
	}

	template <E1_SERIES_TP_DECL>
	POISSON_SERIES sin(const POISSON_SERIES &p)
	{
		POISSON_SERIES retval = p.sin();
		return retval;
	}

	template <E1_SERIES_TP_DECL>
	POISSON_SERIES pow(const POISSON_SERIES &x, const double &y)
	{
		POISSON_SERIES retval(x.pow(y));
		return retval;
	}

	template <E1_SERIES_TP_DECL>
	POISSON_SERIES pow(const POISSON_SERIES &x, const piranha::max_fast_int &n)
	{
		POISSON_SERIES retval(x.pow(n));
		return retval;
	}

	template <E1_SERIES_TP_DECL>
	COMPLEX_POISSON_SERIES pow(const COMPLEX_POISSON_SERIES &x, const double &y)
	{
		COMPLEX_POISSON_SERIES retval(x.pow(y));
		return retval;
	}

	template <E1_SERIES_TP_DECL>
	COMPLEX_POISSON_SERIES pow(const COMPLEX_POISSON_SERIES &x, const piranha::max_fast_int &n)
	{
		COMPLEX_POISSON_SERIES retval(x.pow(n));
		return retval;
	}
}

#endif
