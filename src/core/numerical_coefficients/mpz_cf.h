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

#include <boost/type_traits/integral_constant.hpp> // For lightweight attribute.
#include <cmath>
#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../settings.h"
#include "../type_traits.h" // For lightweight attribute.

namespace piranha
{
	/// Mpz numerical coefficient.
	/**
	 * Arbitrary-size integer coefficient type, to be used as coefficient in piranha::base_series.
	 */
	class mpz_cf: public numerical_container<mpz_class, mpz_cf>
	{
			// Alias for the parent class.
			typedef numerical_container<mpz_class, mpz_cf> ancestor;
		public:
			typedef mpz_class numerical_type;
			// Ctors
			NUMERICAL_CONTAINER_CTORS(mpz_cf);
			// Override norm and evaluation.
			template <class ArgsTuple>
			double norm_(const ArgsTuple &) const {
				return std::abs(m_value.get_d());
			}
			template <class ArgsTuple>
			double eval_(const double &, const ArgsTuple &) const {
				return m_value.get_d();
			}
			// Override swapping for increased efficiency.
			void swap(mpz_cf &other) {
				mpz_swap(m_value.get_mpz_t(),other.m_value.get_mpz_t());
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
			// Multiply and add.
			template <class ArgsTuple>
			void addmul(const mpz_cf &x1, const mpz_cf &x2, const ArgsTuple &) {
				mpz_addmul(m_value.get_mpz_t(), x1.m_value.get_mpz_t(), x2.m_value.get_mpz_t());
			}
			// To shield against GMP exceptions, we need to filter out some cases here.
			template <class ArgsTuple>
			mpz_cf pow_(const int &n, const ArgsTuple &) const {
				mpz_cf retval;
				// If negative, only 1^-something is reasonable.
				if (n < 0) {
					if (m_value == 0) {
						throw division_by_zero();
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						throw unsuitable("Cannot raise integer coefficient different from unity to "
							"negative integer power.");
					}
				} else {
					mpz_pow_ui(retval.m_value.get_mpz_t(), m_value.get_mpz_t(), (size_t)n);
				}
				return retval;
			}
			template <class ArgsTuple>
			mpz_cf pow_(const double &y, const ArgsTuple &) const {
				mpz_cf retval;
				// If negative, only 1^-something is reasonable.
				if (y < 0) {
					if (m_value == 0) {
						throw division_by_zero();
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						throw unsuitable("Cannot raise integer coefficient different from unity to negative real power.");
					}
					// If y == 0, then x**0 == 1 for every x.
				} else if (y == 0) {
					retval.m_value = 1;
					// If y > 0, we can accept only 0^y and 1^y.
				} else {
					if (m_value == 0) {
						retval.m_value = 0;
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						throw unsuitable("Cannot raise integer coefficient different from unity to positive real power.");
					}
				}
				return retval;
			}
			template <class ArgsTuple>
			mpz_cf root_(const int &n_, const ArgsTuple &) const {
				mpz_cf retval;
				if (n_ == 0) {
					throw division_by_zero();
				} else if (n_ == 1) {
					retval = *this;
					return retval;
				} else if (n_ < 0) {
					throw unsuitable("Integer coefficients different from unity cannot be arguments of negative root.");
				}
				const size_t n = static_cast<size_t>(n_);
				if (!mpz_root(retval.m_value.get_mpz_t(),m_value.get_mpz_t(),n)) {
					throw unsuitable("Integer coefficient is not an exact nth root.");
				}
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
			using complex_toolbox::divide_by;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			typedef piranha::mpz_cf value_type;
			NUMERICAL_CONTAINER_CTORS(complex);
			COMPLEX_NUMERICAL_CONTAINER_CTORS;
			// Override numerical container's division to catch divide by zero.
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
			// Override swapping for increased efficiency.
			void swap(complex &other) {
				mpz_swap(m_value.real().get_mpz_t(),other.m_value.real().get_mpz_t());
				mpz_swap(m_value.imag().get_mpz_t(),other.m_value.imag().get_mpz_t());
			}
			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (m_value.real() == 0 && m_value.imag() == 0);
			}
			template <class ArgsTuple>
			double norm_(const ArgsTuple &) const {
				return std::abs(complex<double>(m_value.real().get_d(), m_value.imag().get_d()));
			}
			template <class ArgsTuple>
			complex<double> eval_(const double &, const ArgsTuple &) const {
				return complex<double>(m_value.real().get_d(), m_value.imag().get_d());
			}
			template <class ArgsTuple>
			complex pow_(const int &n, const ArgsTuple &) const {
				complex retval;
				// For negative powers, we must guard against division by zero.
				if (n < 0) {
					if (m_value.real() == 0 && m_value.imag() == 0) {
						throw piranha::division_by_zero();
					} else if (m_value.real() == 1 && m_value.imag() == 0) {
						retval.m_value.real() = 1;
						retval.m_value.imag() = 0;
					} else {
						throw piranha::unsuitable("Cannot raise complex integer coefficient different from unity to "
							"negative integer power.");
					}
				} else {
					retval.m_value.real() = 1;
					retval.m_value.imag() = 0;
					const size_t count = static_cast<size_t>(n);
					for (size_t i = 0; i < count; ++i) {
						retval.m_value *= m_value;
					}
				}
				return retval;
			}
			template <class ArgsTuple>
			complex pow_(const double &y, const ArgsTuple &) const {
				complex retval;
				if (y < 0) {
					if (m_value.real() == 0 && m_value.imag() == 0) {
						throw piranha::division_by_zero();
					} else if (m_value.real() == 1 && m_value.imag() == 0) {
						retval.m_value.real() = 1;
						retval.m_value.imag() = 0;
					} else {
						throw piranha::unsuitable("Cannot raise complex integer coefficient different from unity to "
							"negative real power.");
					}
					// If y == 0, then x**0 == 1 for every x.
				} else if (y == 0) {
					retval.m_value.real() = 1;
					retval.m_value.imag() = 0;
					// If y > 0, we can accept only 0^y and 1^y.
				} else {
					if (m_value.real() == 0 && m_value.imag() == 0) {
						retval.m_value.real() = 0;
						retval.m_value.imag() = 0;
					} else if (m_value.real() == 1 && m_value.imag() == 0) {
						retval.m_value.real() = 1;
						retval.m_value.imag() = 0;
					} else {
						throw piranha::unsuitable("Cannot raise complex integer coefficient different from unity to "
							"positive real power.");
					}
				}
				return retval;
			}
	};
}

namespace piranha
{
	// Specialise lightweight type trait for mpz coefficients.
	template <>
	struct is_lightweight<mpz_cf>: boost::true_type {};

	template <>
	struct is_lightweight<std::complex<mpz_cf> >: boost::true_type {};
}

#endif
