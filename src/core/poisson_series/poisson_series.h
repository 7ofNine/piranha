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
#include "../fourier_series/fourier_series_term.h"
#include "../harmonic_series/base_harmonic_series.h"
#include "../harmonic_series/named_harmonic_series.h"
#include "../mp.h"
#include "../poisson_series_common/common_poisson_series_toolbox.h"
#include "../poisson_series_common/celmec_toolbox.h"
#include "../poisson_series_common/jacobi_anger_toolbox.h"
#include "../polynomial_cf/polynomial_cf.h"

#define POISSON_SERIES                   E1_SERIES(piranha::poisson_series)
#define POISSON_SERIES_POLYNOMIAL_CF     E1_SERIES_COEFFICIENT(piranha::PolynomialCf)
#define POISSON_SERIES_TERM              E1_SERIES_TERM(piranha::FourierSeriesTerm, POISSON_SERIES_POLYNOMIAL_CF)
#define POISSON_SERIES_BASE_ANCESTOR     E1_SERIES_BASE_ANCESTOR(piranha::FourierSeriesTerm, POISSON_SERIES_POLYNOMIAL_CF, POISSON_SERIES)
#define POISSON_SERIES_NAMED_ANCESTOR    E1_SERIES_NAMED_ANCESTOR(piranha::PolyArgsDescr, piranha::TrigArgsDescriptor, POISSON_SERIES_TERM, POISSON_SERIES)
#define POISSON_SERIES_BINOMIAL_ANCESTOR piranha::BinomialExponentiation< POISSON_SERIES>
#define POISSON_SERIES_DEGREE            typename POISSON_SERIES_TERM::CfType::TermType::KeyType::DegreeType
#define POISSON_SERIES_H_DEGREE          typename POISSON_SERIES_TERM::KeyType::HarmonicDegreeType

namespace piranha
{
	template < E1_SERIES_TP_DECL>
	class poisson_series:
				public POISSON_SERIES_BASE_ANCESTOR,
				public POISSON_SERIES_NAMED_ANCESTOR,
				public POISSON_SERIES_BINOMIAL_ANCESTOR,
				public BaseHarmonicSeries<1, 1, POISSON_SERIES_H_DEGREE, POISSON_SERIES>,
				public NamedHarmonicSeries<POISSON_SERIES_H_DEGREE, POISSON_SERIES>,
				public SeriesMultiplication< POISSON_SERIES, Mult1, Trunc1>,
				public JacobiAnger<1, POISSON_SERIES>,
				public CommonPoissonSeries< POISSON_SERIES>,
				public BasePowerSeries<0, 0, POISSON_SERIES_DEGREE, POISSON_SERIES>,
				public NamedPowerSeries< POISSON_SERIES_DEGREE, POISSON_SERIES>,
				public BaseSeriesSpecialFunctions< POISSON_SERIES>,
				public NamedSeriesSpecialFunctions< POISSON_SERIES>,
				public celmec< POISSON_SERIES>,
				boost::ring_operators < POISSON_SERIES,
				boost::ring_operators < POISSON_SERIES, double,
				boost::dividable      < POISSON_SERIES, double,
				boost::ring_operators < POISSON_SERIES, mp_rational,
				boost::dividable      < POISSON_SERIES, mp_rational,
				boost::ring_operators < POISSON_SERIES, mp_integer,
				boost::dividable      < POISSON_SERIES, mp_integer
				> > > > > > >
	{
		public:
			using POISSON_SERIES_BINOMIAL_ANCESTOR::realPower;
			using POISSON_SERIES_BINOMIAL_ANCESTOR::negativeIntegerPower;
			using POISSON_SERIES_BINOMIAL_ANCESTOR::rationalPower;
			using CommonPoissonSeries< POISSON_SERIES >::sub;
			NAMED_SERIES_BOILERPLATE(poisson_series, 0);
	};
}

#define COMPLEX_POISSON_SERIES std::complex<POISSON_SERIES>

#define COMPLEX_POISSON_SERIES_POLYNOMIAL_CF piranha::PolynomialCf<Cf, Key0, Mult0, Trunc0>

#define COMPLEX_POISSON_SERIES_TERM piranha::FourierSeriesTerm<std::complex<COMPLEX_POISSON_SERIES_POLYNOMIAL_CF>, Key1, '|'>

#define COMPLEX_POISSON_SERIES_BASE_ANCESTOR piranha::BaseSeries<COMPLEX_POISSON_SERIES_TERM,'\n', COMPLEX_POISSON_SERIES>

#define COMPLEX_POISSON_SERIES_NAMED_ANCESTOR piranha::NamedSeries<boost::tuple<piranha::PolyArgsDescr, piranha::TrigArgsDescriptor>, \
	COMPLEX_POISSON_SERIES_TERM, COMPLEX_POISSON_SERIES>

#define COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX piranha::BaseSeriesComplex< POISSON_SERIES>

#define COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex< POISSON_SERIES>

#define COMPLEX_POISSON_SERIES_BINOMIAL_ANCESTOR piranha::BinomialExponentiation< COMPLEX_POISSON_SERIES>

#define COMPLEX_POISSON_SERIES_DEGREE typename COMPLEX_POISSON_SERIES_TERM::CfType::TermType::KeyType::DegreeType

#define COMPLEX_POISSON_SERIES_H_DEGREE typename COMPLEX_POISSON_SERIES_TERM::KeyType::HarmonicDegreeType

namespace std
{
	template <E1_SERIES_TP_DECL>
	class complex<POISSON_SERIES>:
				public COMPLEX_POISSON_SERIES_BASE_ANCESTOR,
				public COMPLEX_POISSON_SERIES_NAMED_ANCESTOR,
				public COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX,
				public COMPLEX_POISSON_SERIES_BINOMIAL_ANCESTOR,
				public piranha::BaseHarmonicSeries<1, 1, COMPLEX_POISSON_SERIES_H_DEGREE, COMPLEX_POISSON_SERIES>,
				public piranha::NamedHarmonicSeries<COMPLEX_POISSON_SERIES_H_DEGREE, COMPLEX_POISSON_SERIES>,
				public piranha::SeriesMultiplication< COMPLEX_POISSON_SERIES, Mult1, Trunc1>,
				public piranha::CommonPoissonSeries< COMPLEX_POISSON_SERIES>,
				public piranha::BasePowerSeries<0, 0, COMPLEX_POISSON_SERIES_DEGREE, COMPLEX_POISSON_SERIES>,
				public piranha::NamedPowerSeries< COMPLEX_POISSON_SERIES_DEGREE, COMPLEX_POISSON_SERIES>,
				public piranha::BaseSeriesSpecialFunctions< COMPLEX_POISSON_SERIES>,
				public piranha::NamedSeriesSpecialFunctions< COMPLEX_POISSON_SERIES>,
				boost::ring_operators < COMPLEX_POISSON_SERIES,
				boost::ring_operators < COMPLEX_POISSON_SERIES, double,
				boost::dividable      < COMPLEX_POISSON_SERIES, double,
				boost::ring_operators < COMPLEX_POISSON_SERIES, piranha::mp_rational,
				boost::dividable      < COMPLEX_POISSON_SERIES, piranha::mp_rational,
				boost::ring_operators < COMPLEX_POISSON_SERIES, piranha::mp_integer,
				boost::dividable      < COMPLEX_POISSON_SERIES, piranha::mp_integer,
				boost::ring_operators < COMPLEX_POISSON_SERIES, POISSON_SERIES,
				boost::ring_operators < COMPLEX_POISSON_SERIES, complex<double>,
				boost::dividable      < COMPLEX_POISSON_SERIES, complex<double>
				> > > > > > > > > >
	{
		public:
			using COMPLEX_POISSON_SERIES_BINOMIAL_ANCESTOR::realPower;
			using COMPLEX_POISSON_SERIES_BINOMIAL_ANCESTOR::negativeIntegerPower;
			using COMPLEX_POISSON_SERIES_BINOMIAL_ANCESTOR::rationalPower;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::baseInvert;
			using piranha::CommonPoissonSeries< COMPLEX_POISSON_SERIES>::sub;
			// Ctors.
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(POISSON_SERIES);
	};
}

#endif
