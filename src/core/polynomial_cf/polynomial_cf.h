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
#include <vector>

#include "../base_classes/base_series.h"
#include "../base_classes/base_series_complex_toolbox.h"
#include "../base_classes/base_series_special_functions.h"
#include "../base_classes/cf_series.h"
#include "../base_classes/power_series.h"
#include "../base_classes/series_multiplication.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../none.h"
#include "../polynomial_common/monomial.h"
#include "common_polynomial_cf_toolbox.h"

#define POLYNOMIAL_CF_TERM CF_SERIES_TERM(piranha::monomial,'!')
#define POLYNOMIAL_CF E0_SERIES(piranha::polynomial_cf)
#define POLYNOMIAL_CF_BASE_ANCESTOR CF_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial_cf,'!',',')
#define POLYNOMIAL_CF_CF_ANCESTOR piranha::cf_series< POLYNOMIAL_CF >
#define POLYNOMIAL_CF_MULT_ANCESTOR piranha::series_multiplication< POLYNOMIAL_CF, Multiplier, Truncator>
#define POLYNOMIAL_CF_POWER_SERIES_ANCESTOR power_series<0, 1, POLYNOMIAL_CF >
#define POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR piranha::base_series_special_functions< POLYNOMIAL_CF >
#define POLYNOMIAL_CF_COMMON_ANCESTOR piranha::common_polynomial_cf_toolbox< POLYNOMIAL_CF >

namespace piranha
{
	// NOTE: this assumes that exponents are in position 0 of arguments tuple.
	template <E0_SERIES_TP_DECL>
	class polynomial_cf:
				public POLYNOMIAL_CF_BASE_ANCESTOR,
				public POLYNOMIAL_CF_CF_ANCESTOR,
				public POLYNOMIAL_CF_COMMON_ANCESTOR,
				public POLYNOMIAL_CF_POWER_SERIES_ANCESTOR,
				public POLYNOMIAL_CF_MULT_ANCESTOR,
				public POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR
	{
			typedef POLYNOMIAL_CF_CF_ANCESTOR cf_ancestor;
			typedef POLYNOMIAL_CF_BASE_ANCESTOR base_ancestor;
			typedef POLYNOMIAL_CF_COMMON_ANCESTOR common_ancestor;
			friend class POLYNOMIAL_CF_CF_ANCESTOR;
			friend class POLYNOMIAL_CF_BASE_ANCESTOR;
			friend class POLYNOMIAL_CF_MULT_ANCESTOR;
			// Specify we will use the real_power from the common polynomial cf toolbox.
			using POLYNOMIAL_CF_COMMON_ANCESTOR::real_power;
			using POLYNOMIAL_CF_COMMON_ANCESTOR::negative_integer_power;
			using POLYNOMIAL_CF_COMMON_ANCESTOR::nth_root;
		public:
			using POLYNOMIAL_CF_BASE_ANCESTOR::mult_by;
			using POLYNOMIAL_CF_CF_ANCESTOR::mult_by;
			// Needed typedefs.
			typedef typename Multiplier::template get_type<polynomial_cf, polynomial_cf, none, Truncator> multiplier_type;
			CF_SERIES_CTORS(polynomial_cf);
			template <class ArgsTuple>
			explicit polynomial_cf(const psym_p &p, const int &n, const ArgsTuple &a) {
				base_ancestor::construct_from_psym_p(p, n, a);
			}
	};
}

#define COMPLEX_POLYNOMIAL_CF_TERM COMPLEX_CF_SERIES_TERM(piranha::monomial,'!')
#define COMPLEX_POLYNOMIAL_CF COMPLEX_E0_SERIES(piranha::polynomial_cf)
#define COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR COMPLEX_CF_SERIES_BASE_ANCESTOR(piranha::monomial,piranha::polynomial_cf,'!',',')
#define COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR piranha::cf_series< COMPLEX_POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_MULT_ANCESTOR piranha::series_multiplication< COMPLEX_POLYNOMIAL_CF, Multiplier, Truncator>
#define COMPLEX_POLYNOMIAL_CF_POWER_SERIES_ANCESTOR piranha::power_series<0, 1, COMPLEX_POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX piranha::base_series_complex_toolbox< POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR piranha::base_series_special_functions< COMPLEX_POLYNOMIAL_CF >
#define COMPLEX_POLYNOMIAL_CF_COMMON_ANCESTOR piranha::common_polynomial_cf_toolbox< COMPLEX_POLYNOMIAL_CF >

namespace std
{
	template < E0_SERIES_TP_DECL >
	class complex<POLYNOMIAL_CF>:
				public COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_COMMON_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_POWER_SERIES_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_MULT_ANCESTOR,
				public COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX,
				public COMPLEX_POLYNOMIAL_CF_SPECIAL_FUNCTION_ANCESTOR
	{
			typedef COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR cf_ancestor;
			typedef COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR base_ancestor;
			typedef COMPLEX_POLYNOMIAL_CF_COMMON_ANCESTOR common_ancestor;
			typedef COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX base_complex_toolbox;
			friend class COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR;
			friend class COMPLEX_POLYNOMIAL_CF_MULT_ANCESTOR;
			friend class piranha::base_series_complex_toolbox<POLYNOMIAL_CF>;
			// Specify we will use the real_power from the common polynomial cf toolbox.
			using COMPLEX_POLYNOMIAL_CF_COMMON_ANCESTOR::real_power;
			using COMPLEX_POLYNOMIAL_CF_COMMON_ANCESTOR::negative_integer_power;
			using COMPLEX_POLYNOMIAL_CF_COMMON_ANCESTOR::nth_root;
		public:
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::add;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::add;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::subtract;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::subtract;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::mult_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::mult_by;
			using COMPLEX_POLYNOMIAL_CF_CF_ANCESTOR::mult_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX::divide_by;
			using COMPLEX_POLYNOMIAL_CF_BASE_ANCESTOR::divide_by;
			// Needed typedefs.
			typedef typename Multiplier::template get_type<complex, complex, piranha::none, Truncator> multiplier_type;
			CF_SERIES_CTORS(complex);
			COMPLEX_CF_SERIES_CTORS(COMPLEX_POLYNOMIAL_CF_BASE_COMPLEX_TOOLBOX);
			template <class ArgsTuple>
			explicit complex(const piranha::psym_p &p, const int &n, const ArgsTuple &a) {
				base_ancestor::construct_from_psym_p(p, n, a);
			}
	};
}

#endif
