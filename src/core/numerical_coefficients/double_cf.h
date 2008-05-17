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
	 * This class can be used as coefficient in Poisson series. It encapsulate a double precision
	 * numerical value and provides the means to access it.
	 *
	 * A set of operators is provided to enable interoperability with basic numerical data types.
	 */
	class double_cf: public numerical_container<double, double_cf>
	{
			// Alias for the parent class.
			typedef numerical_container<double, double_cf> ancestor;
		public:
			// Ctors.
			NUMERICAL_CONTAINER_CTORS(double_cf);
			max_fast_int get_int() const {
				const max_fast_int retval((max_fast_int)nearbyint(ancestor::m_value));
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
				//retval.s_value()=math::besselJ(n,g_value());
				//return retval;
				return 0.;
			}
			template <class ArgsTuple>
			double_cf pow(const double &y, const ArgsTuple &) const {
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
		public:
			typedef piranha::double_cf value_type;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			using ancestor::divide_by;
			using complex_toolbox::divide_by;
// Check if these are needed.
			/*      using ancestor::swap;
			      using ancestor::print_plain;
			      using ancestor::print_latex;
			      using ancestor::checkup;
			      using ancestor::invert_sign;
			      using ancestor::t_eval;
			      using ancestor::norm;
			      using ancestor::is_ignorable;
			      using ancestor::is_insertable;
			      using ancestor::needs_padding;
			      using ancestor::pad_right;
			      using ancestor::apply_layout;
			      using ancestor::add;
			      using ancestor::subtract;
			      using ancestor::mult_by;
			      using ancestor::divide_by;
			      using complex_toolbox::mult_by;
			      using complex_toolbox::divide_by;
			      using complex_toolbox::real;
			      using complex_toolbox::imag;
			      using complex_toolbox::set_real;
			      using complex_toolbox::set_imag;*/
			// Start implementation of basic pseries coefficient interface.
			//------------
			NUMERICAL_CONTAINER_CTORS(complex);
			COMPLEX_NUMERICAL_CONTAINER_CTORS;
			// End implementation of complex basic pseries coefficient interface.
			//------------
			template <class ArgsTuple>
			complex pow(const double &y, const ArgsTuple &) const {
				complex retval;
				retval.m_value = std::pow(ancestor::m_value, y);
				return retval;
			}
	};
}
#endif
