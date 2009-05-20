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

#ifndef PIRANHA_DOUBLE_CF_H
#define PIRANHA_DOUBLE_CF_H

#include <boost/numeric/conversion/converter.hpp>
#include <boost/type_traits/integral_constant.hpp> // For lightweight attribute.
#include <cmath>
#include <complex>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../math.h"
#include "../psym.h"
#include "../settings.h" // Numerical zero.
#include "../type_traits.h" // For lightweight attribute.
#include "../utils.h"

namespace piranha
{
	/// Double-precision numerical coefficient.
	/**
	 * This class is meant to be used as coefficient in series. It encapsulate a double precision
	 * numerical value and provides the means for manipulation in the context of a series.
	 */
	class double_cf: public numerical_container<double, double_cf>
	{
			// Alias for the parent class.
			typedef numerical_container<double, double_cf> ancestor;
		public:
			// Ctors.
			NUMERICAL_CONTAINER_CTORS(double_cf,.to_double())
			int get_int() const {
				if (!utils::is_integer(m_value)) {
					piranha_throw(value_error,"cannot convert double coefficient to integer");
				}
				return (int)m_value;
			}
			/// Bessel function of the first kind.
			/**
			 * Uses Boost's math toolkit.
			 */
			template <class ArgsTuple>
			double_cf besselJ(const int &n, const ArgsTuple &) const {
				double_cf retval;
				retval.m_value = piranha::besselJ(n, m_value);
				return retval;
			}
			template <class ArgsTuple>
			std::complex<double_cf> ei(const ArgsTuple &) const;
	};
}

namespace std
{
	template <>
	class complex<piranha::double_cf>:
				public piranha::numerical_container<std::complex<double>, complex<piranha::double_cf> >,
				public piranha::numerical_container_complex_toolbox<piranha::double_cf>
	{
			typedef piranha::numerical_container<std::complex<double>, complex<piranha::double_cf> > ancestor;
			typedef piranha::numerical_container_complex_toolbox<piranha::double_cf> complex_toolbox;
			friend class piranha::numerical_container_complex_toolbox<piranha::double_cf>;
			friend class piranha::double_cf;
		public:
			typedef piranha::double_cf value_type;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			using ancestor::divide_by;
			using complex_toolbox::divide_by;
			NUMERICAL_CONTAINER_CTORS(complex,.to_double())
			COMPLEX_NUMERICAL_CONTAINER_CTORS(.to_complex_double())
			template <class ArgsTuple>
			complex besselJ(const int &n, const ArgsTuple &) {
				// COMPILER_BUG: same as above :/
				complex retval;
				std::complex<double> tmp(piranha::besselJ(n,m_value));
				retval.m_value = tmp;
				return retval;
			}
	};
}

namespace piranha
{
	// Place it here, since double_cf is not a template class and hence this definition
	// relies on the complex specialization for double_cf being already available.
	template <class ArgsTuple>
	inline std::complex<double_cf> double_cf::ei(const ArgsTuple &) const
	{
		std::complex<double_cf> retval;
		retval.m_value = std::polar(1., m_value);
		return retval;
	}

	// Specialise lightweight type trait for double coefficients.
	template <>
	struct is_lightweight<double_cf>: boost::true_type {};

	template <>
	struct is_lightweight<std::complex<double_cf> >: boost::true_type {};
}

#endif
