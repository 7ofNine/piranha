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
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../integer_typedefs.h"
#include "../polynomial_common/monomial.h"
#include "../polynomial_common/common_polynomial_toolbox.h"

#define POLYNOMIAL_TERM E0_SERIES_TERM(piranha::monomial)
#define POLYNOMIAL E0_SERIES(piranha::polynomial)
#define POLYNOMIAL_BASE_ANCESTOR E0_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial)
#define POLYNOMIAL_NAMED_ANCESTOR E0_SERIES_NAMED_ANCESTOR(boost::tuple<poly_args_descr>,piranha::polynomial)
#define POLYNOMIAL_MULT_ANCESTOR piranha::series_multiplication< POLYNOMIAL, Multiplier, Truncator>
#define POLYNOMIAL_COMMON_ANCESTOR piranha::common_polynomial_toolbox< POLYNOMIAL >
#define POLYNOMIAL_POWER_SERIES_ANCESTOR power_series<0,1,POLYNOMIAL >
#define POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR piranha::base_series_special_functions< POLYNOMIAL >
#define POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR piranha::named_series_special_functions< POLYNOMIAL >

namespace piranha
{
	template < E0_SERIES_TP_DECL = std::allocator<char> >
	class polynomial:
				public POLYNOMIAL_BASE_ANCESTOR,
				public POLYNOMIAL_NAMED_ANCESTOR,
				public POLYNOMIAL_POWER_SERIES_ANCESTOR,
				public POLYNOMIAL_MULT_ANCESTOR,
				public POLYNOMIAL_COMMON_ANCESTOR,
				public POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR,
				public POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR,
				boost::ring_operators < POLYNOMIAL,
				boost::ring_operators < POLYNOMIAL, max_fast_int,
				boost::ring_operators < POLYNOMIAL, double,
				boost::dividable < POLYNOMIAL, max_fast_int,
				boost::dividable < POLYNOMIAL, double
				> > > > >
	{
			typedef POLYNOMIAL_NAMED_ANCESTOR named_ancestor;
			typedef POLYNOMIAL_BASE_ANCESTOR base_ancestor;
			friend class POLYNOMIAL_NAMED_ANCESTOR;
			friend class POLYNOMIAL_BASE_ANCESTOR;
			friend class POLYNOMIAL_MULT_ANCESTOR;
			friend class POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class named_series_complex_toolbox<POLYNOMIAL>;
			// Override power functions with the ones from the common polynomial toolbox.
			using POLYNOMIAL_COMMON_ANCESTOR::real_power;
			using POLYNOMIAL_COMMON_ANCESTOR::negative_integer_power;
			using POLYNOMIAL_COMMON_ANCESTOR::nth_root;
			using POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR::besselJ;
			using POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR::dbesselJ;
		public:
			using named_ancestor::norm;
			using POLYNOMIAL_COMMON_ANCESTOR::norm;
			using named_ancestor::pow;
			using base_ancestor::pow;
			using named_ancestor::root;
			using base_ancestor::root;
			using named_ancestor::partial;
			using base_ancestor::partial;
			using POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR::besselJ;
			using POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR::dbesselJ;
			// This is needed because called from the norm override in poly toolbox.
			using base_ancestor::eval;
			// Needed typedefs.
			typedef typename Multiplier::template get_type < polynomial, polynomial,
			typename named_ancestor::args_tuple_type, Truncator > multiplier_type;
			// Boilerplate.
			NAMED_SERIES_BOILERPLATE(polynomial, 0);
	};
}

#define COMPLEX_POLYNOMIAL_TERM COMPLEX_E0_SERIES_TERM(piranha::monomial)
#define COMPLEX_POLYNOMIAL COMPLEX_E0_SERIES(piranha::polynomial)
#define COMPLEX_POLYNOMIAL_BASE_ANCESTOR COMPLEX_E0_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial)
#define COMPLEX_POLYNOMIAL_NAMED_ANCESTOR COMPLEX_E0_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::poly_args_descr>, \
		piranha::polynomial)
#define COMPLEX_POLYNOMIAL_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_POLYNOMIAL, Multiplier, Truncator>
#define COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX piranha::base_series_complex_toolbox<POLYNOMIAL>
#define COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex_toolbox<POLYNOMIAL>
#define COMPLEX_POLYNOMIAL_COMMON_ANCESTOR piranha::common_polynomial_toolbox< COMPLEX_POLYNOMIAL >
#define COMPLEX_POLYNOMIAL_POWER_SERIES_ANCESTOR piranha::power_series<0,1,COMPLEX_POLYNOMIAL >
#define COMPLEX_POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR piranha::base_series_special_functions< COMPLEX_POLYNOMIAL >
#define COMPLEX_POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR piranha::named_series_special_functions< COMPLEX_POLYNOMIAL >

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL>:
				public COMPLEX_POLYNOMIAL_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_ANCESTOR,
				public COMPLEX_POLYNOMIAL_MULT_ANCESTOR,
				public COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_COMMON_ANCESTOR,
				public COMPLEX_POLYNOMIAL_POWER_SERIES_ANCESTOR,
				public COMPLEX_POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR,
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
			typedef COMPLEX_POLYNOMIAL_NAMED_ANCESTOR named_ancestor;
			typedef COMPLEX_POLYNOMIAL_BASE_ANCESTOR base_ancestor;
			typedef COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX base_complex_toolbox;
			friend class COMPLEX_POLYNOMIAL_NAMED_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_BASE_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_MULT_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX;
			friend class COMPLEX_POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX;
			// Override power_functions with the ones from the common polynomial toolbox.
			using COMPLEX_POLYNOMIAL_COMMON_ANCESTOR::real_power;
			using COMPLEX_POLYNOMIAL_COMMON_ANCESTOR::negative_integer_power;
			using COMPLEX_POLYNOMIAL_COMMON_ANCESTOR::nth_root;
			using COMPLEX_POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR::besselJ;
			using COMPLEX_POLYNOMIAL_BASE_SPECIAL_FUNCTIONS_ANCESTOR::dbesselJ;
		public:
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::real;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::real;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::imag;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::imag;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::construct_from_real;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::construct_from_real;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::construct_from_real_imag;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::construct_from_real_imag;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::add;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::add;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::subtract;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::mult_by;
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_POLYNOMIAL_BASE_ANCESTOR::divide_by;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator+=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator-=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator*=;
			using COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_POLYNOMIAL_NAMED_ANCESTOR::operator/=;
			using named_ancestor::norm;
			using COMPLEX_POLYNOMIAL_COMMON_ANCESTOR::norm;
			using named_ancestor::pow;
			using base_ancestor::pow;
			using named_ancestor::root;
			using base_ancestor::root;
			using named_ancestor::partial;
			using base_ancestor::partial;
			using base_ancestor::eval;
			using COMPLEX_POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR::besselJ;
			using COMPLEX_POLYNOMIAL_NAMED_SPECIAL_FUNCTIONS_ANCESTOR::dbesselJ;
			// Needed typedefs.
			typedef typename Multiplier::template get_type < complex, complex,
			typename named_ancestor::args_tuple_type, Truncator > multiplier_type;
			// Boilerplate and additional ctors.
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX);
	};
}

// Overload standard math functions for polynomials.
namespace std
{
	template < E0_SERIES_TP_DECL >
	POLYNOMIAL pow(const POLYNOMIAL &x, const piranha::max_fast_int &n)
	{
		POLYNOMIAL retval(x.pow(n));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	POLYNOMIAL pow(const POLYNOMIAL &x, const double &y)
	{
		POLYNOMIAL retval(x.pow(y));
		return retval;
	}

	template < E0_SERIES_TP_DECL >
	COMPLEX_POLYNOMIAL pow(const COMPLEX_POLYNOMIAL &x, const piranha::max_fast_int &n)
	{
		COMPLEX_POLYNOMIAL retval(x.pow(n));
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
