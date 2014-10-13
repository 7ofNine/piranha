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

#ifndef PIRANHA_MPQ_CF_H
#define PIRANHA_MPQ_CF_H

#include <boost/type_traits/integral_constant.hpp>
#include <cmath>
#include <complex>

#include "../base_classes/numerical_container.h"
#include "../base_classes/numerical_container_complex_toolbox.h"
#include "../mp.h"
#include "../type_traits.h"

namespace piranha
{
	/// Rational numerical coefficient for series.
	/**
	 * Arbitrary-size rational coefficient type, to be used as coefficient in piranha::BaseSeries.
	 */
	class mpq_cf: public NumericalContainer<mp_rational, mpq_cf>
	{
			// Alias for the parent class.
			typedef NumericalContainer<mp_rational, mpq_cf> ancestor;

		public:

			explicit mpq_cf(): ancestor() {}

			template <class T, class ArgsTuple>
			explicit mpq_cf(const T &x, const ArgsTuple &argsTuple): ancestor(x, argsTuple) {}

			template <class ArgsTuple>
			explicit mpq_cf(const Psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}

			/// Override print in Tex mode.
			template <class ArgsTuple>
			void print_tex(std::ostream &outStream, const ArgsTuple &) const
			{
				get_value().print_tex(outStream);
			}

			// Override norm and evaluation.
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

			int to_int() const 
            {
				return get_value().to_int();
			}

			template <class ArgsTuple>
			std::complex<mpq_cf> ei(const ArgsTuple &) const;
	};

	/// is_ring_exact type trait specialisation for mpq_cf.
	template <>
	struct is_ring_exact<mpq_cf>: boost::true_type {};

	/// is_divint_exact type trait specialisation for mpq_cf.
	template <>
	struct is_divint_exact<mpq_cf>: boost::true_type {};
}


namespace std
{
	template <>
	class complex<piranha::mpq_cf>:	public piranha::NumericalContainer<complex<piranha::mp_rational>, complex<piranha::mpq_cf> >,
				                    public piranha::NumericalContainerComplexToolbox<piranha::mpq_cf>
	{
			typedef piranha::NumericalContainer<complex<piranha::mp_rational>, complex<piranha::mpq_cf> > ancestor;
			typedef piranha::NumericalContainerComplexToolbox<piranha::mpq_cf> complex_toolbox;

		public:

			explicit complex(): ancestor() {}

			template <class T, class ArgsTuple>
			explicit complex(const T &x, const ArgsTuple &argsTuple): ancestor(x, argsTuple) {}

			template <class ArgsTuple>
			explicit complex(const piranha::Psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}

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

			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (get_value() == 0);
			}
	};
}


namespace piranha
{
	template <class ArgsTuple>
	inline std::complex<mpq_cf> mpq_cf::ei(const ArgsTuple &) const
	{
		PIRANHA_THROW(value_error,"rational coefficient is unsuitable for complex exponentiation");
	}
}

#endif
