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

#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include <boost/operators.hpp>
#include <boost/type_traits/integral_constant.hpp>

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/base_power_series.h"
#include "../base_classes/named_power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../mp.h"
#include "../polynomial_common/base_polynomial.h"
#include "../polynomial_common/monomial.h"
#include "../type_traits.h"
#include "named_polynomial.h"

#define POLYNOMIAL_TERM              E0_SERIES_TERM(piranha::Monomial)
#define POLYNOMIAL                   E0_SERIES(piranha::Polynomial)
#define POLYNOMIAL_BASE_ANCESTOR     E0_SERIES_BASE_ANCESTOR(piranha::Monomial, piranha::Polynomial)
#define POLYNOMIAL_NAMED_ANCESTOR    E0_SERIES_NAMED_ANCESTOR(boost::tuple<PolyArgsDescr>, POLYNOMIAL_TERM, piranha::Polynomial)
#define POLYNOMIAL_BINOMIAL_ANCESTOR piranha::BinomialExponentiation<POLYNOMIAL>
#define POLYNOMIAL_DEGREE            typename POLYNOMIAL_TERM::KeyType::DegreeType
#define POLYNOMIAL_BASE_POLYNOMIAL_ANCESTOR  piranha::BasePolynomial<0, POLYNOMIAL >
#define POLYNOMIAL_NAMED_POLYNOMIAL_ANCESTOR piranha::NamedPolynomial<POLYNOMIAL >

namespace piranha
{
	/// Polynomial class.
	template < E0_SERIES_TP_DECL = std::allocator<char> >
	class Polynomial:
				public POLYNOMIAL_BASE_ANCESTOR,     //piranha::BaseSeries<piranha::Monomial<Cf,Key,'|',Allocator>),'\n', Allocator, piranha::Polynomial<Cf,Key,Multiplier,Truncator,Allocator> >
				public POLYNOMIAL_NAMED_ANCESTOR,
				public POLYNOMIAL_BINOMIAL_ANCESTOR, //piranha::binomial_exponentiation<piranha::Polynomial<Cf,Key,Multiplier,Truncator,Allocator> >
				public POLYNOMIAL_BASE_POLYNOMIAL_ANCESTOR, //piranha::BasePolynomial<0, piranha::Polynomial<Cf,Key,Multiplier,Truncator,Allocator> >
				public POLYNOMIAL_NAMED_POLYNOMIAL_ANCESTOR, //piranha::NamedPolynomial<piranha::Polynomial<Cf,Key,Multiplier,Truncator,Allocator> >
				public BasePowerSeries<0, 1, POLYNOMIAL_DEGREE, POLYNOMIAL>,
				public NamedPowerSeries<POLYNOMIAL_DEGREE, POLYNOMIAL>,
				public series_multiplication< POLYNOMIAL, Multiplier, Truncator>,
				public BaseSeriesSpecialFunctions< POLYNOMIAL>,
				public named_series_special_functions< POLYNOMIAL>,
				boost::ring_operators < POLYNOMIAL,
				boost::ring_operators < POLYNOMIAL, double,
				boost::dividable      < POLYNOMIAL, double,
				boost::ring_operators < POLYNOMIAL, mp_rational,
				boost::dividable      < POLYNOMIAL, mp_rational,
				boost::ring_operators < POLYNOMIAL, mp_integer,
				boost::dividable      < POLYNOMIAL, mp_integer
				> > > > > > >
	{
		public:

			using POLYNOMIAL_BINOMIAL_ANCESTOR::realPower;
			using POLYNOMIAL_BINOMIAL_ANCESTOR::negativeIntegerPower;
			using POLYNOMIAL_BINOMIAL_ANCESTOR::rationalPower;
			// Boilerplate.
			NAMED_SERIES_BOILERPLATE(Polynomial, 0);
	};
}

#define COMPLEX_POLYNOMIAL_TERM COMPLEX_E0_SERIES_TERM(piranha::Monomial)
#define COMPLEX_POLYNOMIAL COMPLEX_E0_SERIES(piranha::Polynomial)
#define COMPLEX_POLYNOMIAL_BASE_ANCESTOR COMPLEX_E0_SERIES_BASE_ANCESTOR(piranha::Monomial, piranha::Polynomial)
#define COMPLEX_POLYNOMIAL_NAMED_ANCESTOR COMPLEX_E0_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::PolyArgsDescr>, \
		COMPLEX_POLYNOMIAL_TERM, piranha::Polynomial)
#define COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX piranha::BaseSeriesComplex<POLYNOMIAL>
#define COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex<POLYNOMIAL>
#define COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR piranha::BinomialExponentiation< COMPLEX_POLYNOMIAL>
#define COMPLEX_POLYNOMIAL_DEGREE typename COMPLEX_POLYNOMIAL_TERM::KeyType::DegreeType
#define COMPLEX_POLYNOMIAL_BASE_POLYNOMIAL_ANCESTOR piranha::BasePolynomial<0, COMPLEX_POLYNOMIAL>
#define COMPLEX_POLYNOMIAL_NAMED_POLYNOMIAL_ANCESTOR piranha::NamedPolynomial<COMPLEX_POLYNOMIAL>

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL>:
				public COMPLEX_POLYNOMIAL_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_ANCESTOR,
				public COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_NAMED_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR,
				public COMPLEX_POLYNOMIAL_BASE_POLYNOMIAL_ANCESTOR,
				public COMPLEX_POLYNOMIAL_NAMED_POLYNOMIAL_ANCESTOR,
				public piranha::series_multiplication< COMPLEX_POLYNOMIAL, Multiplier, Truncator>,
				public piranha::BasePowerSeries<0, 1, COMPLEX_POLYNOMIAL_DEGREE, COMPLEX_POLYNOMIAL>,
				public piranha::NamedPowerSeries< COMPLEX_POLYNOMIAL_DEGREE, COMPLEX_POLYNOMIAL>,
				public piranha::BaseSeriesSpecialFunctions< COMPLEX_POLYNOMIAL>,
				public piranha::named_series_special_functions< COMPLEX_POLYNOMIAL>,
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
		public:
			using COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR::realPower;
			using COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR::negativeIntegerPower;
			using COMPLEX_POLYNOMIAL_BINOMIAL_ANCESTOR::rationalPower;
			// Complex overrides.
			using COMPLEX_POLYNOMIAL_BASE_COMPLEX_TOOLBOX::baseInvert;
			// Boilerplate and additional ctors.
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(POLYNOMIAL);
	};
}

#endif
