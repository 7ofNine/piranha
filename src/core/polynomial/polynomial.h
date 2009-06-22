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
#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/common_comparisons.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/base_power_series.h"
#include "../base_classes/named_power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/toolbox.h"
#include "../mp.h"
#include "../polynomial_common/monomial.h"
#include "../polynomial_common/base_polynomial_toolbox.h"

#define POLYNOMIAL_TERM E0_SERIES_TERM(piranha::monomial)
#define POLYNOMIAL E0_SERIES(piranha::polynomial)
#define POLYNOMIAL_BASE_ANCESTOR E0_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial)
#define POLYNOMIAL_NAMED_ANCESTOR E0_SERIES_NAMED_ANCESTOR(boost::tuple<poly_args_descr>, POLYNOMIAL_TERM, piranha::polynomial)
#define POLYNOMIAL_BINOMIAL_ANCESTOR piranha::toolbox<piranha::binomial_exponentiation< POLYNOMIAL ,piranha::term_key_degree_comparison> >
#define POLYNOMIAL_DEGREE typename POLYNOMIAL_TERM::key_type::degree_type

namespace piranha
{
	template < E0_SERIES_TP_DECL = std::allocator<char> >
	class polynomial:
				public POLYNOMIAL_BASE_ANCESTOR,
				public POLYNOMIAL_NAMED_ANCESTOR,
				public POLYNOMIAL_BINOMIAL_ANCESTOR,
				public toolbox<base_power_series<0,1,POLYNOMIAL_DEGREE,POLYNOMIAL > >,
				public toolbox<named_power_series< POLYNOMIAL_DEGREE,POLYNOMIAL > >,
				public toolbox<series_multiplication< POLYNOMIAL, Multiplier, Truncator> >,
				public toolbox<base_polynomial_toolbox< POLYNOMIAL > >,
				public toolbox<base_series_special_functions< POLYNOMIAL > >,
				public toolbox<named_series_special_functions< POLYNOMIAL > >,
				boost::ring_operators < POLYNOMIAL,
				boost::ring_operators < POLYNOMIAL, double,
				boost::dividable < POLYNOMIAL, double,
				boost::ring_operators < POLYNOMIAL, mp_rational,
				boost::dividable < POLYNOMIAL, mp_rational,
				boost::ring_operators < POLYNOMIAL, mp_integer,
				boost::dividable < POLYNOMIAL, mp_integer
				> > > > > > >
	{
			template <class>
			friend class toolbox;
			using POLYNOMIAL_BINOMIAL_ANCESTOR::real_power;
			using POLYNOMIAL_BINOMIAL_ANCESTOR::negative_integer_power;
			using POLYNOMIAL_BINOMIAL_ANCESTOR::rational_power;
			using toolbox<base_polynomial_toolbox< POLYNOMIAL > >::base_norm;
		public:
			// Boilerplate.
			NAMED_SERIES_BOILERPLATE(polynomial, 0);
	};
}

#define COMPLEX_POLYNOMIAL_TERM COMPLEX_E0_SERIES_TERM(piranha::monomial)
#define COMPLEX_POLYNOMIAL COMPLEX_E0_SERIES(piranha::polynomial)
#define COMPLEX_POLYNOMIAL_BASE_ANCESTOR COMPLEX_E0_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial)
#define COMPLEX_POLYNOMIAL_NAMED_ANCESTOR COMPLEX_E0_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::poly_args_descr>, \
		COMPLEX_POLYNOMIAL_TERM, piranha::polynomial)
#define COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX piranha::toolbox<piranha::base_series_complex<POLYNOMIAL> >
#define COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX piranha::toolbox<piranha::named_series_complex<POLYNOMIAL> >
#define COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR piranha::toolbox<piranha::binomial_exponentiation< COMPLEX_POLYNOMIAL , piranha::term_key_degree_comparison> >
#define COMPLEX_POLYNOMIAL_DEGREE typename COMPLEX_POLYNOMIAL_TERM::key_type::degree_type

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL>:
				public COMPLEX_POLYNOMIAL_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_ANCESTOR,
				public COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR,
				public piranha::toolbox<piranha::series_multiplication< COMPLEX_POLYNOMIAL, Multiplier, Truncator> >,
				public piranha::toolbox<piranha::base_polynomial_toolbox< COMPLEX_POLYNOMIAL > >,
				public piranha::toolbox<piranha::base_power_series<0,1,COMPLEX_POLYNOMIAL_DEGREE,COMPLEX_POLYNOMIAL > >,
				public piranha::toolbox<piranha::named_power_series< COMPLEX_POLYNOMIAL_DEGREE,COMPLEX_POLYNOMIAL > >,
				public piranha::toolbox<piranha::base_series_special_functions< COMPLEX_POLYNOMIAL > >,
				public piranha::toolbox<piranha::named_series_special_functions< COMPLEX_POLYNOMIAL > >,
				boost::ring_operators < COMPLEX_POLYNOMIAL,
				boost::ring_operators < COMPLEX_POLYNOMIAL, double,
				boost::dividable < COMPLEX_POLYNOMIAL, double,
				boost::ring_operators < COMPLEX_POLYNOMIAL, piranha::mp_rational,
				boost::dividable < COMPLEX_POLYNOMIAL, piranha::mp_rational,
				boost::ring_operators < COMPLEX_POLYNOMIAL, piranha::mp_integer,
				boost::dividable < COMPLEX_POLYNOMIAL, piranha::mp_integer,
				boost::ring_operators < COMPLEX_POLYNOMIAL, POLYNOMIAL,
				boost::ring_operators < COMPLEX_POLYNOMIAL, complex<double>,
				boost::dividable < COMPLEX_POLYNOMIAL, complex<double>
				> > > > > > > > > >
	{
			template <class>
			friend class piranha::toolbox;
			using COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR::real_power;
			using COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR::negative_integer_power;
			using COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR::rational_power;
			// Complex overrides.
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::base_inv;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::base_add;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::base_add;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::base_subtract;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::base_subtract;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::base_mult_by;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::base_mult_by;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::base_divide_by;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::base_divide_by;
			// Polynomial norm override.
			using piranha::toolbox<piranha::base_polynomial_toolbox< COMPLEX_POLYNOMIAL > >::base_norm;
		public:
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator==;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator==;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator!=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator!=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator+=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator-=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator*=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator/=;
			// Boilerplate and additional ctors.
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(POLYNOMIAL);
	};
}

// Overload standard math functions for polynomials.
namespace std
{
	template < E0_SERIES_TP_DECL >
	POLYNOMIAL pow(const POLYNOMIAL &x, const double &y)
	{
		POLYNOMIAL retval(x.pow(y));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	COMPLEX_POLYNOMIAL pow(const COMPLEX_POLYNOMIAL &x, const double &y)
	{
		COMPLEX_POLYNOMIAL retval(x.pow(y));
		return retval;
	}
}

#endif
