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
#include <boost/type_traits/integral_constant.hpp>
#include <cmath>
#include <complex>
#include <memory> // For default allocator.

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../mp.h"
#include "../harmonic_series/base_harmonic_series.h"
#include "../harmonic_series/named_harmonic_series.h"
#include "../poisson_series_common/jacobi_anger_toolbox.h"
#include "base_fourier_series.h"
#include "common_fourier_series_toolbox.h"
#include "fourier_series_term.h"
#include "named_fourier_series.h"

#define FOURIER_SERIES_TERM E0_SERIES_TERM(piranha::FourierSeriesTerm)
#define FOURIER_SERIES E0_SERIES(piranha::fourier_series)
#define FOURIER_SERIES_BASE_ANCESTOR E0_SERIES_BASE_ANCESTOR(piranha::FourierSeriesTerm,piranha::fourier_series)
#define FOURIER_SERIES_NAMED_ANCESTOR E0_SERIES_NAMED_ANCESTOR(boost::tuple<trig_args_descr>, FOURIER_SERIES_TERM ,piranha::fourier_series)
#define FOURIER_SERIES_BINOMIAL_ANCESTOR piranha::binomial_exponentiation< FOURIER_SERIES>
#define FOURIER_SERIES_H_DEGREE typename FOURIER_SERIES_TERM::key_type::h_degree_type
#define FOURIER_SERIES_BASE_FOURIER_SERIES_ANCESTOR piranha::BaseFourierSeries<0,FOURIER_SERIES>
#define FOURIER_SERIES_NAMED_FOURIER_SERIES_ANCESTOR piranha::named_fourier_series<FOURIER_SERIES>

namespace piranha
{
	template < E0_SERIES_TP_DECL = std::allocator<char> >
	class fourier_series:
				public FOURIER_SERIES_BASE_ANCESTOR,
				public FOURIER_SERIES_NAMED_ANCESTOR,
				public FOURIER_SERIES_BINOMIAL_ANCESTOR,
				public FOURIER_SERIES_BASE_FOURIER_SERIES_ANCESTOR,
				public FOURIER_SERIES_NAMED_FOURIER_SERIES_ANCESTOR,
				public BaseHarmonicSeries<0,1,FOURIER_SERIES_H_DEGREE,FOURIER_SERIES>,
				public named_harmonic_series<FOURIER_SERIES_H_DEGREE,FOURIER_SERIES>,
				public common_fourier_series< FOURIER_SERIES>,
				public series_multiplication< FOURIER_SERIES, Multiplier, Truncator>,
				public base_series_special_functions< FOURIER_SERIES>,
				public named_series_special_functions< FOURIER_SERIES>,
				public jacobi_anger<0, FOURIER_SERIES>,
				boost::ring_operators < FOURIER_SERIES,
				boost::ring_operators < FOURIER_SERIES, double,
				boost::dividable < FOURIER_SERIES, double,
				boost::ring_operators < FOURIER_SERIES, mp_rational,
				boost::dividable < FOURIER_SERIES, mp_rational,
				boost::ring_operators < FOURIER_SERIES, mp_integer,
				boost::dividable < FOURIER_SERIES, mp_integer
				> > > > > > >
	{
		public:
			using FOURIER_SERIES_BINOMIAL_ANCESTOR::real_power;
			using FOURIER_SERIES_BINOMIAL_ANCESTOR::negative_integer_power;
			using FOURIER_SERIES_BINOMIAL_ANCESTOR::rational_power;
			// Boilerplate
			NAMED_SERIES_BOILERPLATE(fourier_series, 0);
	};
}

#define COMPLEX_FOURIER_SERIES_TERM COMPLEX_E0_SERIES_TERM(piranha::FourierSeriesTerm)
#define COMPLEX_FOURIER_SERIES COMPLEX_E0_SERIES(piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_BASE_ANCESTOR COMPLEX_E0_SERIES_BASE_ANCESTOR(piranha::FourierSeriesTerm,piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR COMPLEX_E0_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::trig_args_descr>, \
		COMPLEX_FOURIER_SERIES_TERM , piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX piranha::BaseSeriesComplex< FOURIER_SERIES>
#define COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX piranha::named_series_complex< FOURIER_SERIES>
#define COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR piranha::binomial_exponentiation< COMPLEX_FOURIER_SERIES>
#define COMPLEX_FOURIER_SERIES_H_DEGREE typename COMPLEX_FOURIER_SERIES_TERM::key_type::h_degree_type
#define COMPLEX_FOURIER_SERIES_BASE_FOURIER_SERIES_ANCESTOR piranha::BaseFourierSeries<0,COMPLEX_FOURIER_SERIES>
#define COMPLEX_FOURIER_SERIES_NAMED_FOURIER_SERIES_ANCESTOR piranha::named_fourier_series<COMPLEX_FOURIER_SERIES>

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<FOURIER_SERIES>:
				public COMPLEX_FOURIER_SERIES_BASE_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX,
				public COMPLEX_FOURIER_SERIES_BASE_FOURIER_SERIES_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_NAMED_FOURIER_SERIES_ANCESTOR,
				public piranha::BaseHarmonicSeries<0,1,COMPLEX_FOURIER_SERIES_H_DEGREE,COMPLEX_FOURIER_SERIES>,
				public piranha::named_harmonic_series<COMPLEX_FOURIER_SERIES_H_DEGREE,COMPLEX_FOURIER_SERIES>,
				public piranha::series_multiplication< COMPLEX_FOURIER_SERIES, Multiplier, Truncator>,
				public piranha::common_fourier_series < COMPLEX_FOURIER_SERIES>,
				public piranha::base_series_special_functions< COMPLEX_FOURIER_SERIES>,
				public piranha::named_series_special_functions< COMPLEX_FOURIER_SERIES>,
				public COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR,
				boost::ring_operators < COMPLEX_FOURIER_SERIES,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, double,
				boost::dividable < COMPLEX_FOURIER_SERIES, double,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, piranha::mp_rational,
				boost::dividable < COMPLEX_FOURIER_SERIES, piranha::mp_rational,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, piranha::mp_integer,
				boost::dividable < COMPLEX_FOURIER_SERIES, piranha::mp_integer,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, FOURIER_SERIES,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, complex<double>,
				boost::dividable < COMPLEX_FOURIER_SERIES, complex<double>
				> > > > > > > > > >
	{
		public:
			using COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR::real_power;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR::negative_integer_power;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR::rational_power;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::base_inv;
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(FOURIER_SERIES);
	};
}

#endif
