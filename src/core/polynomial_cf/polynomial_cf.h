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
#include "../base_classes/toolbox.h"
#include "../exceptions.h"
#include "../polynomial_common/monomial.h"
#include "common_polynomial_cf_toolbox.h"

#define POLYNOMIAL_CF_TERM CF_SERIES_TERM(piranha::monomial,'!')
#define POLYNOMIAL_CF E0_SERIES(piranha::polynomial_cf)
#define POLYNOMIAL_CF_BASE_ANCESTOR CF_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial_cf,'!','?')
#define POLYNOMIAL_CF_CF_ANCESTOR piranha::toolbox<piranha::cf_series< POLYNOMIAL_CF_TERM, POLYNOMIAL_CF > >
#define POLYNOMIAL_CF_BINOMIAL_ANCESTOR piranha::toolbox<piranha::binomial_exponentiation< POLYNOMIAL_CF > >
#define POLYNOMIAL_CF_DEGREE typename POLYNOMIAL_CF_TERM::key_type::degree_type

namespace piranha
{
	// NOTE: this assumes that exponents are in position 0 of arguments tuple. Can be made configurable at a later stage.
	template <E0_SERIES_TP_DECL>
	class polynomial_cf:
				public POLYNOMIAL_CF_BASE_ANCESTOR,
				public POLYNOMIAL_CF_CF_ANCESTOR,
				public POLYNOMIAL_CF_BINOMIAL_ANCESTOR,
				public toolbox<common_polynomial_cf< POLYNOMIAL_CF > >,
				public toolbox<base_power_series<0, 1, POLYNOMIAL_CF_DEGREE, POLYNOMIAL_CF > >,
				public toolbox<cf_power_series< POLYNOMIAL_CF_DEGREE, POLYNOMIAL_CF > >,
				public toolbox<series_multiplication< POLYNOMIAL_CF, Multiplier, Truncator> >,
				public toolbox<base_series_special_functions< POLYNOMIAL_CF > >,
				public toolbox<cf_series_special_functions< POLYNOMIAL_CF > >
	{
			template <class>
			friend class toolbox;
			// Specify we will use the power functions from the binomial toolbox.
			using POLYNOMIAL_CF_BINOMIAL_ANCESTOR::real_power;
			using POLYNOMIAL_CF_BINOMIAL_ANCESTOR::negative_integer_power;
			using POLYNOMIAL_CF_BINOMIAL_ANCESTOR::rational_power;
		public:
			CF_SERIES_CTORS(polynomial_cf);
			template <class ArgsTuple>
			explicit polynomial_cf(const psym &p, const int &n, const ArgsTuple &a) {
				this->base_construct_from_psym(p, n, a);
			}
	};
}

#define COMPLEX_POLYNOMIAL_CF_TERM COMPLEX_CF_SERIES_TERM(piranha::monomial,'!')
#define COMPLEX_POLYNOMIAL_CF COMPLEX_E0_SERIES(piranha::polynomial_cf)
#define COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR COMPLEX_CF_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial_cf,'!','?')
#define COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR piranha::toolbox<piranha::cf_series< COMPLEX_POLYNOMIAL_CF_TERM, COMPLEX_POLYNOMIAL_CF > >
#define COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX piranha::toolbox<piranha::base_series_complex< POLYNOMIAL_CF > >
#define COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX piranha::toolbox<piranha::cf_series_complex< POLYNOMIAL_CF > >
#define COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR piranha::toolbox<piranha::binomial_exponentiation< COMPLEX_POLYNOMIAL_CF > >
#define COMPLEX_POLYNOMIAL_CF_DEGREE typename COMPLEX_POLYNOMIAL_CF_TERM::key_type::degree_type

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL_CF>:
				public COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR,
				public piranha::toolbox<piranha::common_polynomial_cf< COMPLEX_POLYNOMIAL_CF > >,
				public piranha::toolbox<piranha::base_power_series<0, 1, COMPLEX_POLYNOMIAL_CF_DEGREE, COMPLEX_POLYNOMIAL_CF > >,
				public piranha::toolbox<piranha::cf_power_series< COMPLEX_POLYNOMIAL_CF_DEGREE, COMPLEX_POLYNOMIAL_CF > >,
				public piranha::toolbox<piranha::series_multiplication< COMPLEX_POLYNOMIAL_CF, Multiplier, Truncator> >,
				public piranha::toolbox<piranha::base_series_special_functions< COMPLEX_POLYNOMIAL_CF > >,
				public piranha::toolbox<piranha::cf_series_special_functions< COMPLEX_POLYNOMIAL_CF > >
	{
			template <class>
			friend class piranha::toolbox;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR::real_power;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR::negative_integer_power;
			using COMPLEX_POLYNOMIAL_CF_BINOMIAL_ANCESTOR::rational_power;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::base_inv;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::base_add;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::base_add;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::base_subtract;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::base_subtract;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::base_mult_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::base_mult_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::base_divide_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::base_divide_by;
		public:
			using COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX::operator==;
			using COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR::operator==;
			using COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR::mult_by;
			using COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR::divide_by;
			using COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX::add;
			using COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR::add;
			using COMPLEX_POLYNOMIAL_CF_CF_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR::subtract;
			CF_SERIES_CTORS(complex);
			COMPLEX_CF_SERIES_CTORS(POLYNOMIAL_CF);
			template <class ArgsTuple>
			explicit complex(const piranha::psym &p, const int &n, const ArgsTuple &a) {
				this->base_construct_from_psym(p, n, a);
			}
	};
}

#endif
