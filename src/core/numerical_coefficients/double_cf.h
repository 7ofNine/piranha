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
			// Ctors and dtor.
			/// Empty constructor.
			explicit double_cf(): ancestor::numerical_container() {}
			/// Constructor from string.
			template <class ArgsTuple>
			explicit double_cf(const std::string &s, const ArgsTuple &a): ancestor::numerical_container(s, a) {}
			/// Constructor from integer.
			template <class ArgsTuple>
			explicit double_cf(const max_fast_int &val, const ArgsTuple &a): ancestor::numerical_container(val, a) {}
			/// Constructor from double.
			template <class ArgsTuple>
			explicit double_cf(const double &val, const ArgsTuple &a): ancestor::numerical_container(val, a) {}
			/// Constructor from psym.
			template <class ArgsTuple>
			explicit double_cf(const psym_p &p, const int &n, const ArgsTuple &a): ancestor::numerical_container(p, n, a) {}
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
			typedef complex self;
			friend class piranha::numerical_container_complex_toolbox<piranha::double_cf>;
		public:
			typedef piranha::double_cf value_type;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
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
			// Basic ctors and dtor.
			explicit complex(): ancestor::numerical_container() {}
			template <class ArgsTuple>
			explicit complex(const std::string &s, const ArgsTuple &a): ancestor::numerical_container(s, a) {}
			template <class ArgsTuple>
			explicit complex(const piranha::max_fast_int &n, const ArgsTuple &a): ancestor::numerical_container(n, a) {}
			template <class ArgsTuple>
			explicit complex(const double &x, const ArgsTuple &a): ancestor::numerical_container(x, a) {}
			complex(const complex &c): ancestor::numerical_container(c), complex_toolbox::numerical_container_complex_toolbox(c) {}
			// Complex specific contructors.
			template <class ArgsTuple>
			explicit complex(const std::complex<piranha::max_fast_int> &c, const ArgsTuple &):
				complex_toolbox::numerical_container_complex_toolbox(c) {}
			template <class ArgsTuple>
			explicit complex(const std::complex<double> &c, const ArgsTuple &):
				complex_toolbox::numerical_container_complex_toolbox(c) {}
			template <class ArgsTuple>
			explicit complex(const value_type &r, const ArgsTuple &): complex_toolbox::numerical_container_complex_toolbox(r) {}
			template <class ArgsTuple>
			explicit complex(const value_type &r, const value_type &i, const ArgsTuple &):
					complex_toolbox::numerical_container_complex_toolbox(r, i) {}
			// Operators.
			complex &operator=(const self &val2) {
				return assign_self(val2);
			}
			complex &operator=(const value_type &r2) {
				return assign_self(r2);
			}
			// End implementation of complex basic pseries coefficient interface.
			//------------
	};
}
#endif
