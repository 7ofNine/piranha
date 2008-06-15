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
#include <cmath>
#include <complex>

#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h" // besselJ.
#include "../psym.h"
#include "../settings.h" // Numerical zero.
#include "../base_classes/numerical_container.h"

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
			typedef double numerical_type;
			// Ctors.
			NUMERICAL_CONTAINER_CTORS(double_cf);
			max_fast_int get_int() const {
				typedef boost::numeric::converter<max_fast_int,double> double_to_int;
				const max_fast_int retval((max_fast_int)double_to_int::nearbyint(ancestor::m_value));
				if (std::abs(ancestor::m_value - retval) > settings::numerical_zero()) {
					throw(unsuitable("Cannot convert double coefficient to integer."));
				}
				return retval;
			}
			/// Bessel function of the first kind.
			/**
			 * Uses Boost's math toolkit.
			 */
			template <class ArgsTuple>
			double_cf besselJ(int n, const ArgsTuple &) const {
				double_cf retval;
				retval.m_value = piranha::besselJ(n, m_value);
				return retval;
			}
			template <class ArgsTuple>
			double_cf pow(const max_fast_int &n, const ArgsTuple &) const {
				return pow_helper(n);
			}
			template <class ArgsTuple>
			double_cf pow(const double &y, const ArgsTuple &) const {
				return pow_helper(y);
			}
			template <class ArgsTuple>
			std::complex<double_cf> complexp(const ArgsTuple &) const;
		private:
			template <class Number>
			double_cf pow_helper(const Number &y) const {
				double_cf retval;
				retval.m_value = std::pow(ancestor::m_value, y);
				return retval;
			}
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
			NUMERICAL_CONTAINER_CTORS(complex);
			COMPLEX_NUMERICAL_CONTAINER_CTORS;
			template <class ArgsTuple>
			complex pow(const piranha::max_fast_int &n, const ArgsTuple &) const {
				return pow_helper(n);
			}
			template <class ArgsTuple>
			complex pow(const double &y, const ArgsTuple &) const {
				return pow_helper(y);
			}
		private:
			template <class Number>
			complex pow_helper(const Number &y) const {
				complex retval;
				retval.m_value = std::pow(ancestor::m_value, y);
				return retval;
			}
	};
}

namespace piranha
{
	// Place it here, since double_cf is not a template class and hence this definition
	// relies on the complex specialization for double_cf being already available.
	template <class ArgsTuple>
	inline std::complex<double_cf> double_cf::complexp(const ArgsTuple &) const
	{
		std::complex<double_cf> retval;
		retval.m_value = std::polar(1., m_value);
		return retval;
	}
}

#endif
