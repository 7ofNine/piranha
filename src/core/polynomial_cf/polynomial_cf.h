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

#include <boost/type_traits/integral_constant.hpp>
#include <complex>

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/cf_series.h"
#include "../base_classes/cf_series_complex_toolbox.h"
#include "../base_classes/cf_series_special_functions.h"
#include "../base_classes/base_power_series.h"
#include "../base_classes/cf_power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../exceptions.h"
#include "../polynomial_common/base_polynomial.h"
#include "../polynomial_common/monomial.h"
#include "../type_traits.h"
#include "common_polynomial_cf_toolbox.h"

#define POLYNOMIAL_CF_TERM CF_SERIES_TERM(piranha::Monomial,'!')
#define POLYNOMIAL_CF E0_SERIES(piranha::PolynomialCf)
#define POLYNOMIAL_CF_BASE_ANCESTOR CF_SERIES_BASE_ANCESTOR(piranha::Monomial, piranha::PolynomialCf, '!', '?')
#define POLYNOMIAL_CF_CF_ANCESTOR piranha::CfSeries< POLYNOMIAL_CF_TERM, POLYNOMIAL_CF>
#define POLYNOMIAL_CF_BINOMIAL_ANCESTOR piranha::BinomialExponentiation< POLYNOMIAL_CF>
#define POLYNOMIAL_CF_DEGREE typename POLYNOMIAL_CF_TERM::KeyType::DegreeType
#define POLYNOMIAL_CF_BASE_POLYNOMIAL_ANCESTOR piranha::BasePolynomial<0, POLYNOMIAL_CF>

namespace piranha
{
	// NOTE: this assumes that exponents are in position 0 of arguments tuple. Can be made configurable at a later stage.
	template <E0_SERIES_TP_DECL>
	class PolynomialCf:
				public POLYNOMIAL_CF_BASE_ANCESTOR,
				public POLYNOMIAL_CF_CF_ANCESTOR,
				public POLYNOMIAL_CF_BINOMIAL_ANCESTOR,
				public POLYNOMIAL_CF_BASE_POLYNOMIAL_ANCESTOR,
				public CommonPolynomialCf< POLYNOMIAL_CF>,
				public BasePowerSeries<0, 1, POLYNOMIAL_CF_DEGREE, POLYNOMIAL_CF>,
				public CfPowerSeries< POLYNOMIAL_CF_DEGREE, POLYNOMIAL_CF>,
				public SeriesMultiplication< POLYNOMIAL_CF, Multiplier, Truncator>,
				public BaseSeriesSpecialFunctions< POLYNOMIAL_CF>,
				public CfSeriesSpecialFunctions< POLYNOMIAL_CF>
	{
		public:
			// Specify we will use the power functions from the binomial toolbox.
			using POLYNOMIAL_CF_BINOMIAL_ANCESTOR::realPower;
			using POLYNOMIAL_CF_BINOMIAL_ANCESTOR::negativeIntegerPower;
			using POLYNOMIAL_CF_BINOMIAL_ANCESTOR::rationalPower;

			CF_SERIES_CTORS(PolynomialCf);

			template <class ArgsTuple>
			explicit PolynomialCf(const Psym &p, const int &n, const ArgsTuple &a)
            {
				this->baseConstructFromPsym(p, n, a);
			}
	};
}

#define COMPLEX_POLYNOMIAL_CF_TERM COMPLEX_CF_SERIES_TERM(piranha::Monomial, '!')
#define COMPLEX_POLYNOMIAL_CF COMPLEX_E0_SERIES(piranha::PolynomialCf)
#define COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR COMPLEX_CF_SERIES_BASE_ANCESTOR(piranha::Monomial, piranha::PolynomialCf, '!', '?')
#define COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR piranha::CfSeries< COMPLEX_POLYNOMIAL_CF_TERM, COMPLEX_POLYNOMIAL_CF>
#define COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX piranha::BaseSeriesComplex< POLYNOMIAL_CF>
#define COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX piranha::CfSeriesComplex< POLYNOMIAL_CF>
#define COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR piranha::BinomialExponentiation< COMPLEX_POLYNOMIAL_CF>
#define COMPLEX_POLYNOMIAL_CF_DEGREE typename COMPLEX_POLYNOMIAL_CF_TERM::KeyType::DegreeType
#define COMPLEX_POLYNOMIAL_CF_BASE_POLYNOMIAL_ANCESTOR piranha::BasePolynomial<0,COMPLEX_POLYNOMIAL_CF>

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL_CF>:
				public COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_BASE_POLYNOMIAL_ANCESTOR,
				public piranha::CommonPolynomialCf< COMPLEX_POLYNOMIAL_CF>,
				public piranha::BasePowerSeries<0, 1, COMPLEX_POLYNOMIAL_CF_DEGREE, COMPLEX_POLYNOMIAL_CF>,
				public piranha::CfPowerSeries< COMPLEX_POLYNOMIAL_CF_DEGREE, COMPLEX_POLYNOMIAL_CF>,
				public piranha::SeriesMultiplication< COMPLEX_POLYNOMIAL_CF, Multiplier, Truncator>,
				public piranha::BaseSeriesSpecialFunctions< COMPLEX_POLYNOMIAL_CF>,
				public piranha::CfSeriesSpecialFunctions< COMPLEX_POLYNOMIAL_CF>
	{
		public:
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR::realPower;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR::negativeIntegerPower;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR::rationalPower;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::baseInvert;
			CF_SERIES_CTORS(complex);
			COMPLEX_CF_SERIES_CTORS(POLYNOMIAL_CF);
			template <class ArgsTuple>
			explicit complex(const piranha::Psym &p, const int &n, const ArgsTuple &a) {
				this->baseConstructFromPsym(p, n, a);
			}
	};
}

#endif
