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

#ifndef PIRANHA_MPZ_CF_H
#define PIRANHA_MPZ_CF_H

#include <boost/type_traits/integral_constant.hpp>
#include <cmath>
#include <complex>

#include "../base_classes/numerical_container.h"
#include "../base_classes/numerical_container_complex_toolbox.h"
#include "../mp.h"
#include "../type_traits.h"

namespace piranha
{
	/// Multiprecision integer coefficient.
	/**
	 * Arbitrary-size integer coefficient type, to be used as coefficient in piranha::BaseSeries.
	 */
	class mpz_cf: public NumericalContainer<mp_integer, mpz_cf>
	{
			// Alias for the parent class.
			typedef NumericalContainer<mp_integer, mpz_cf> ancestor;

		public:

			explicit mpz_cf(): ancestor() {}

			template <class T, class ArgsTuple>
			explicit mpz_cf(const T &x, const ArgsTuple &args_tuple): ancestor(x,args_tuple) {}

			template <class ArgsTuple>
			explicit mpz_cf(const psym &p, const int &n, const ArgsTuple &a): ancestor(p,n,a) {}

			// Implement norm and evaluation.
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const 
            {
				return std::abs(get_value().to_double());
			}

			template <class ArgsTuple>
			double eval(const double &, const ArgsTuple &) const 
            {
				return get_value().to_double();
			}

			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const 
            {
				return (get_value() == 0);
			}

			template <class ArgsTuple>
			mpz_cf besselJ(const int &n, const ArgsTuple &args_tuple) const
			{
				if (get_value() != 0) 
                {
					piranha_throw(value_error,"cannot compute Bessel function of non-zero integer coefficient");
				}

				return (n == 0) ? mpz_cf(1,args_tuple) : mpz_cf(0,args_tuple);
			}
	};

	/// is_ring_exact type trait specialisation for mpz_cf.
	template <>
	struct is_ring_exact<mpz_cf>: boost::true_type {};
}

namespace std
{
	template <>
	class complex<piranha::mpz_cf>: public piranha::NumericalContainer<complex<piranha::mp_integer>, complex<piranha::mpz_cf> >,
		                            public piranha::NumericalContainerComplexToolbox<piranha::mpz_cf>
	{
			typedef piranha::NumericalContainer<complex<piranha::mp_integer>, complex<piranha::mpz_cf> > ancestor;
			typedef piranha::NumericalContainerComplexToolbox<piranha::mpz_cf> complex_toolbox;

		public:
			explicit complex(): ancestor() {}

			template <class T, class ArgsTuple>
			explicit complex(const T &x, const ArgsTuple &args_tuple): ancestor(x,args_tuple) {}

			template <class ArgsTuple>
			explicit complex(const piranha::psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}

			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const 
            {
				return get_value() == 0;
			}

			template <class ArgsTuple>
			double norm(const ArgsTuple &) const 
            {
				return abs(get_value().to_complex_double());
			}

			template <class ArgsTuple>
			complex<double> eval(const double &, const ArgsTuple &) const 
            {
				return get_value().to_complex_double();
			}
	};
}

#endif
