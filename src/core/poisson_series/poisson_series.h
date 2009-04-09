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
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/named_series.h"
#include "../base_classes/named_series_complex_toolbox.h"
#include "../base_classes/named_series_special_functions.h"
#include "../base_classes/power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../base_classes/toolbox.h"
#include "../fourier_series/fourier_series_term.h"
#include "../poisson_series_common/common_poisson_series_toolbox.h"
#include "../poisson_series_common/celmec_toolbox.h"
#include "../polynomial_cf/polynomial_cf.h"

#define POISSON_SERIES E1_SERIES(piranha::poisson_series)
#define POISSON_SERIES_POLYNOMIAL_CF E1_SERIES_COEFFICIENT(piranha::polynomial_cf)
#define POISSON_SERIES_TERM E1_SERIES_TERM(piranha::fourier_series_term,POISSON_SERIES_POLYNOMIAL_CF)
#define POISSON_SERIES_BASE_ANCESTOR E1_SERIES_BASE_ANCESTOR(piranha::fourier_series_term,POISSON_SERIES_POLYNOMIAL_CF,POISSON_SERIES)
#define POISSON_SERIES_NAMED_ANCESTOR E1_SERIES_NAMED_ANCESTOR(piranha::poly_args_descr, piranha::trig_args_descr, POISSON_SERIES_TERM, POISSON_SERIES)
#define POISSON_SERIES_MULT_ANCESTOR piranha::toolbox<piranha::series_multiplication< POISSON_SERIES, Mult1, Trunc1> >
#define POISSON_SERIES_COMMON_ANCESTOR piranha::toolbox<piranha::common_poisson_series< POISSON_SERIES > >
#define POISSON_SERIES_POWER_SERIES_ANCESTOR piranha::toolbox<piranha::power_series<0, 0, POISSON_SERIES > >
#define POISSON_SERIES_BASE_SPECIAL_FUNCTIONS_ANCESTOR piranha::toolbox<piranha::base_series_special_functions< POISSON_SERIES > >
#define POISSON_SERIES_NAMED_SPECIAL_FUNCTIONS_ANCESTOR piranha::toolbox<piranha::named_series_special_functions< POISSON_SERIES > >
#define POISSON_SERIES_CELMEC_ANCESTOR piranha::toolbox<celmec< POISSON_SERIES > >

namespace piranha
{
	template < E1_SERIES_TP_DECL = std::allocator<char> >
	class poisson_series:
				public POISSON_SERIES_BASE_ANCESTOR,
				public POISSON_SERIES_NAMED_ANCESTOR,
				public POISSON_SERIES_MULT_ANCESTOR,
				public POISSON_SERIES_COMMON_ANCESTOR,
				public POISSON_SERIES_POWER_SERIES_ANCESTOR,
				public POISSON_SERIES_BASE_SPECIAL_FUNCTIONS_ANCESTOR,
				public POISSON_SERIES_NAMED_SPECIAL_FUNCTIONS_ANCESTOR,
				public POISSON_SERIES_CELMEC_ANCESTOR,
				boost::ring_operators < POISSON_SERIES,
				boost::ring_operators < POISSON_SERIES, double,
				boost::dividable < POISSON_SERIES, double
				> > >
	{
			template <class>
			friend class toolbox;
			template <int, int, class>
			friend class trig_array;
			// Additional friendship required by substitution.
			template <int, int, class>
			friend class expo_array;
			using POISSON_SERIES_COMMON_ANCESTOR::real_power;
			using POISSON_SERIES_COMMON_ANCESTOR::negative_integer_power;
			using POISSON_SERIES_COMMON_ANCESTOR::nth_root;
		public:
			using POISSON_SERIES_COMMON_ANCESTOR::sub;
			NAMED_SERIES_BOILERPLATE(poisson_series, 0);
	};
}

#define COMPLEX_POISSON_SERIES std::complex<POISSON_SERIES>
#define COMPLEX_POISSON_SERIES_POLYNOMIAL_CF piranha::polynomial_cf<Cf,Key0,Mult0,Trunc0,Allocator>
#define COMPLEX_POISSON_SERIES_TERM piranha::fourier_series_term<std::complex<COMPLEX_POISSON_SERIES_POLYNOMIAL_CF>,Key1,'|',Allocator>
#define COMPLEX_POISSON_SERIES_BASE_ANCESTOR piranha::toolbox<piranha::base_series<COMPLEX_POISSON_SERIES_TERM,'\n',Allocator,COMPLEX_POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_NAMED_ANCESTOR piranha::toolbox<piranha::named_series<boost::tuple<piranha::poly_args_descr,piranha::trig_args_descr>, \
	COMPLEX_POISSON_SERIES_TERM, COMPLEX_POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_MULT_ANCESTOR piranha::toolbox<piranha::series_multiplication< COMPLEX_POISSON_SERIES, Mult1, Trunc1> >
#define COMPLEX_POISSON_SERIES_COMMON_ANCESTOR piranha::toolbox<piranha::common_poisson_series< COMPLEX_POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_POWER_SERIES_ANCESTOR piranha::toolbox<piranha::power_series<0, 0, COMPLEX_POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX piranha::toolbox<piranha::base_series_complex< POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX piranha::toolbox<piranha::named_series_complex< POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_BASE_SPECIAL_FUNCTIONS_ANCESTOR piranha::toolbox<piranha::base_series_special_functions< COMPLEX_POISSON_SERIES > >
#define COMPLEX_POISSON_SERIES_NAMED_SPECIAL_FUNCTIONS_ANCESTOR piranha::toolbox<piranha::named_series_special_functions< COMPLEX_POISSON_SERIES > >

namespace std
{
	template <E1_SERIES_TP_DECL>
	class complex<POISSON_SERIES>:
				public COMPLEX_POISSON_SERIES_BASE_ANCESTOR,
				public COMPLEX_POISSON_SERIES_NAMED_ANCESTOR,
				public COMPLEX_POISSON_SERIES_MULT_ANCESTOR,
				public COMPLEX_POISSON_SERIES_COMMON_ANCESTOR,
				public COMPLEX_POISSON_SERIES_POWER_SERIES_ANCESTOR,
				public COMPLEX_POISSON_SERIES_BASE_SPECIAL_FUNCTIONS_ANCESTOR,
				public COMPLEX_POISSON_SERIES_NAMED_SPECIAL_FUNCTIONS_ANCESTOR,
				public COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX,
				boost::ring_operators < COMPLEX_POISSON_SERIES,
				boost::ring_operators < COMPLEX_POISSON_SERIES, double,
				boost::dividable < COMPLEX_POISSON_SERIES, double,
				boost::ring_operators < COMPLEX_POISSON_SERIES, POISSON_SERIES,
				boost::ring_operators < COMPLEX_POISSON_SERIES, complex<double>,
				boost::dividable < COMPLEX_POISSON_SERIES, complex<double>
				> > > > > >
	{
			template <class>
			friend class piranha::toolbox;
			using COMPLEX_POISSON_SERIES_COMMON_ANCESTOR::real_power;
			using COMPLEX_POISSON_SERIES_COMMON_ANCESTOR::negative_integer_power;
			using COMPLEX_POISSON_SERIES_COMMON_ANCESTOR::nth_root;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::base_inv;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::base_add;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::base_add;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::base_subtract;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::base_subtract;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::base_mult_by;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::base_mult_by;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::base_divide_by;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::base_divide_by;
		public:
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::operator==;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::operator==;
			using COMPLEX_POISSON_SERIES_BASE_COMPLEX_TOOLBOX::operator!=;
			using COMPLEX_POISSON_SERIES_BASE_ANCESTOR::operator!=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator+=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator+=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator-=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator-=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator*=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator*=;
			using COMPLEX_POISSON_SERIES_NAMED_COMPLEX_TOOLBOX::operator/=;
			using COMPLEX_POISSON_SERIES_NAMED_ANCESTOR::operator/=;
			using COMPLEX_POISSON_SERIES_COMMON_ANCESTOR::sub;
			// Ctors.
			NAMED_SERIES_BOILERPLATE(complex, 0);
			COMPLEX_NAMED_SERIES_CTORS(POISSON_SERIES);
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
	COMPLEX_POISSON_SERIES pow(const COMPLEX_POISSON_SERIES &x, const double &y)
	{
		COMPLEX_POISSON_SERIES retval(x.pow(y));
		return retval;
	}
}

#endif
