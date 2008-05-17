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

#include <cmath>
#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../settings.h"

namespace piranha
{
	/// Mpz numerical coefficient.
	/**
	 * Arbitrary-size integer coefficient type, to be used as coefficient in piranha::base_series.
	 *
	 * A set of operators is provided to enable interoperability with basic numerical data types.
	 */
	class mpz_cf: public numerical_container<mpz_class, mpz_cf>
	{
			// Alias for the parent class.
			typedef numerical_container<mpz_class, mpz_cf> ancestor;
		public:
			// Ctors
			NUMERICAL_CONTAINER_CTORS(mpz_cf);
			// Override norm and evaluation.
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(m_value.get_d());
			}
			// Override division to catch divide by zero.
			template <class ArgsTuple>
			mpz_cf &divide_by(const max_fast_int &n, const ArgsTuple &a) {
				if (n == 0) {
					throw division_by_zero();
				}
				return ancestor::divide_by(n, a);
			}
			template <class ArgsTuple>
			mpz_cf &divide_by(const double &x, const ArgsTuple &a) {
				if (x == 0) {
					throw division_by_zero();
				}
				return ancestor::divide_by(x, a);
			}
			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (m_value == 0);
			}
			template <class ArgsTuple>
			double eval(const double &, const ArgsTuple &) const {
				return m_value.get_d();
			}
			// Multiply and add.
			template <class ArgsTuple>
			void addmul(const mpz_cf &x1, const mpz_cf &x2, const ArgsTuple &) {
				mpz_addmul(m_value.get_mpz_t(), x1.m_value.get_mpz_t(), x2.m_value.get_mpz_t());
			}
			template <class ArgsTuple>
			mpz_cf pow(const double &y, const ArgsTuple &) const {
				mpz_cf retval;
				// If value = 1, then any power is ok, just return 1.
				if (m_value == 1) {
					retval.m_value = 1;
					return retval;
				}
				const max_fast_int pow_n((max_fast_int)nearbyint(y));
				if (std::abs(pow_n - y) > settings::numerical_zero()) {
					throw(unsuitable("Cannot raise integer coefficient different from unity to real power."));
				}
				if (pow_n < 0) {
					throw(unsuitable("Cannot raise integer coefficient different from unity to negative integer power."));
				}
				mpz_pow_ui(retval.m_value.get_mpz_t(), m_value.get_mpz_t(), (size_t)pow_n);
				return retval;
			}
	};
}

namespace std
{
	template <>
	class complex<piranha::mpz_cf>:
				public piranha::numerical_container<std::complex<mpz_class>, complex<piranha::mpz_cf> >,
				public piranha::numerical_container_complex_toolbox<piranha::mpz_cf>
	{
			typedef piranha::numerical_container<std::complex<mpz_class>, complex<piranha::mpz_cf> > ancestor;
			typedef piranha::numerical_container_complex_toolbox<piranha::mpz_cf> complex_toolbox;
			friend class piranha::numerical_container_complex_toolbox<piranha::mpz_cf>;
		public:
			typedef piranha::mpz_cf value_type;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			using ancestor::divide_by;
			using complex_toolbox::divide_by;
			NUMERICAL_CONTAINER_CTORS(complex);
			COMPLEX_NUMERICAL_CONTAINER_CTORS;
			// Override division to catch divide by zero.
			template <class ArgsTuple>
			complex &divide_by(const piranha::max_fast_int &n, const ArgsTuple &a) {
				if (n == 0) {
					throw piranha::division_by_zero();
				}
				return ancestor::divide_by(n, a);
			}
			template <class ArgsTuple>
			complex &divide_by(const double &x, const ArgsTuple &a) {
				if (x == 0) {
					throw piranha::division_by_zero();
				}
				return ancestor::divide_by(x, a);
			}
			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (m_value.real() == 0 and m_value.imag() == 0);
			}
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(complex<double>(m_value.real().get_d(),m_value.imag().get_d()));
			}
			template <class ArgsTuple>
			complex<double> eval(const double &, const ArgsTuple &) const {
				return complex<double>(m_value.real().get_d(),m_value.imag().get_d());
			}
			template <class ArgsTuple>
			complex pow(const double &y, const ArgsTuple &) const {
				complex retval;
				retval.m_value = std::pow(ancestor::m_value, y);
				return retval;
			}
	};
}

#endif
