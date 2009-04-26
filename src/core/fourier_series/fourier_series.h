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
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/binomial_exponentiation_toolbox.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/common_comparisons.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/toolbox.h"
#include "../poisson_series_common/jacobi_anger_toolbox.h"
#include "common_fourier_series_toolbox.h"
#include "fourier_series_term.h"

#define FOURIER_SERIES_TERM E0_SERIES_TERM(piranha::fourier_series_term)
#define FOURIER_SERIES E0_SERIES(piranha::fourier_series)
#define FOURIER_SERIES_BASE_ANCESTOR E0_SERIES_BASE_ANCESTOR(piranha::fourier_series_term,piranha::fourier_series)
#define FOURIER_SERIES_NAMED_ANCESTOR E0_SERIES_NAMED_ANCESTOR(boost::tuple<trig_args_descr>, FOURIER_SERIES_TERM ,piranha::fourier_series)
#define FOURIER_SERIES_BINOMIAL_ANCESTOR piranha::toolbox<piranha::binomial_exponentiation< FOURIER_SERIES, piranha::fs_binomial_sorter > >

namespace piranha
{
	template < E0_SERIES_TP_DECL = std::allocator<char> >
	class fourier_series:
				public FOURIER_SERIES_BASE_ANCESTOR,
				public FOURIER_SERIES_NAMED_ANCESTOR,
				public FOURIER_SERIES_BINOMIAL_ANCESTOR,
				public toolbox<common_fourier_series< FOURIER_SERIES > >,
				public toolbox<series_multiplication< FOURIER_SERIES, Multiplier, Truncator> >,
				public toolbox<base_series_special_functions< FOURIER_SERIES > >,
				public toolbox<named_series_special_functions< FOURIER_SERIES > >,
				public toolbox<jacobi_anger<0, FOURIER_SERIES > >,
				boost::ring_operators < FOURIER_SERIES,
				boost::ring_operators < FOURIER_SERIES, double,
				boost::dividable < FOURIER_SERIES, double
				> > >
	{
			template <class>
			friend class toolbox;
			using FOURIER_SERIES_BINOMIAL_ANCESTOR::real_power;
			using FOURIER_SERIES_BINOMIAL_ANCESTOR::negative_integer_power;
			using FOURIER_SERIES_BINOMIAL_ANCESTOR::nth_root;
		public:
			// Boilerplate
			NAMED_SERIES_BOILERPLATE(fourier_series, 0);
	};
}

#define COMPLEX_FOURIER_SERIES_TERM COMPLEX_E0_SERIES_TERM(piranha::fourier_series_term)
#define COMPLEX_FOURIER_SERIES COMPLEX_E0_SERIES(piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_BASE_ANCESTOR COMPLEX_E0_SERIES_BASE_ANCESTOR(piranha::fourier_series_term,piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR COMPLEX_E0_SERIES_NAMED_ANCESTOR(boost::tuple<piranha::trig_args_descr>, \
		COMPLEX_FOURIER_SERIES_TERM , piranha::fourier_series)
#define COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX piranha::toolbox<piranha::base_series_complex< FOURIER_SERIES > >
#define COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX piranha::toolbox<piranha::named_series_complex< FOURIER_SERIES > >
#define COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR piranha::toolbox<piranha::binomial_exponentiation< COMPLEX_FOURIER_SERIES, piranha::fs_binomial_sorter > >

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<FOURIER_SERIES>:
				public COMPLEX_FOURIER_SERIES_BASE_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR,
				public COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX,
				public piranha::toolbox<piranha::series_multiplication< COMPLEX_FOURIER_SERIES, Multiplier, Truncator> >,
				public piranha::toolbox<piranha::common_fourier_series < COMPLEX_FOURIER_SERIES > >,
				public piranha::toolbox<piranha::base_series_special_functions< COMPLEX_FOURIER_SERIES > >,
				public piranha::toolbox<piranha::named_series_special_functions< COMPLEX_FOURIER_SERIES > >,
				public COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR,
				boost::ring_operators < COMPLEX_FOURIER_SERIES,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, double,
				boost::dividable < COMPLEX_FOURIER_SERIES, double,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, FOURIER_SERIES,
				boost::ring_operators < COMPLEX_FOURIER_SERIES, complex<double>,
				boost::dividable < COMPLEX_FOURIER_SERIES, complex<double>
				> > > > > >
	{
			template <class>
			friend class piranha::toolbox;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR::real_power;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR::negative_integer_power;
			using COMPLEX_FOURIER_SERIES_BINOMIAL_ANCESTOR::nth_root;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::base_inv;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::base_add;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::base_add;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::base_subtract;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::base_subtract;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::base_mult_by;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::base_mult_by;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::base_divide_by;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::base_divide_by;
		public:
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::operator==;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::operator==;
			using COMPLEX_FOURIER_SERIES_BASE_COMPLEX_TOOLBOX::operator!=;
			using COMPLEX_FOURIER_SERIES_BASE_ANCESTOR::operator!=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator+=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator-=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator*=;
			using COMPLEX_FOURIER_SERIES_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_FOURIER_SERIES_NAMED_ANCESTOR::operator/=;
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(FOURIER_SERIES);
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
	COMPLEX_FOURIER_SERIES pow(const COMPLEX_FOURIER_SERIES &x, const double &y)
	{
		COMPLEX_FOURIER_SERIES retval(x.pow(y));
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
