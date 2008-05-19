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

#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../settings.h"

namespace piranha
{
	/// Mpq numerical coefficient.
	/**
	 * Arbitrary-size rational coefficient type, to be used as coefficient in piranha::base_series.
	 *
	 * A set of operators is provided to enable interoperability with basic numerical data types.
	 */
	class mpq_cf: public numerical_container<mpq_class, mpq_cf>
	{
			// Alias for the parent class.
			typedef numerical_container<mpq_class, mpq_cf> ancestor;
		public:
			// Ctors. Do not use macro since we need to canonicalize in one of the ctors.
			/// Empty constructor.
			explicit mpq_cf(): ancestor::numerical_container() {}
			/// Constructor from string.
			template <class ArgsTuple>
			explicit mpq_cf(const std::string &s, const ArgsTuple &a): ancestor::numerical_container(s, a) {
				// We need to canonicalize when reading from string.
				m_value.canonicalize();
			}
			/// Constructor from integer.
			template <class ArgsTuple>
			explicit mpq_cf(const max_fast_int &val, const ArgsTuple &a): ancestor::numerical_container(val, a) {}
			/// Constructor from double.
			template <class ArgsTuple>
			explicit mpq_cf(const double &val, const ArgsTuple &a): ancestor::numerical_container(val, a) {}
			/// Constructor from psym.
			template <class ArgsTuple>
			explicit mpq_cf(const psym_p &p, const int &n, const ArgsTuple &a): ancestor::numerical_container(p, n, a) {}
			// Override norm and evaluation.
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(m_value.get_d());
			}
			template <class ArgsTuple>
			double eval(const double &, const ArgsTuple &) const {
				return m_value.get_d();
			}
			// Override division to catch divide by zero.
			template <class ArgsTuple>
			mpq_cf &divide_by(const max_fast_int &n, const ArgsTuple &a) {
				if (n == 0) {
					throw division_by_zero();
				}
				return ancestor::divide_by(n, a);
			}
			template <class ArgsTuple>
			mpq_cf &divide_by(const double &x, const ArgsTuple &a) {
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
			max_fast_int get_int() const {
				max_fast_int retval = m_value.get_num().get_si();
				if (m_value.get_den() != 1) {
					throw(unsuitable("Cannot convert rational coefficient to integer."));
				}
				return retval;
			}
			template <class ArgsTuple>
			mpq_cf pow(const double &y, const ArgsTuple &) const {
				mpq_cf retval;
				// If value = 1, then any power is ok, just return 1.
				if (m_value == 1) {
					retval.m_value = 1;
				} else {
					const max_fast_int pow_n((max_fast_int)nearbyint(y));
					if (std::abs(pow_n - y) > settings::numerical_zero()) {
						throw(unsuitable("Cannot raise rational coefficient different from unity to real power."));
					}
					if (pow_n < 0) {
						mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (size_t)(-pow_n));
						mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (size_t)(-pow_n));
					} else {
						mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (size_t)pow_n);
						mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (size_t)pow_n);
					}
				}
				return retval;
			}
	};
}

namespace std
{
	template <>
	class complex<piranha::mpq_cf>:
				public piranha::numerical_container<std::complex<mpq_class>, complex<piranha::mpq_cf> >,
				public piranha::numerical_container_complex_toolbox<piranha::mpq_cf>
	{
			typedef piranha::numerical_container<std::complex<mpq_class>, complex<piranha::mpq_cf> > ancestor;
			typedef piranha::numerical_container_complex_toolbox<piranha::mpq_cf> complex_toolbox;
			friend class piranha::numerical_container_complex_toolbox<piranha::mpq_cf>;
		public:
			typedef piranha::mpq_cf value_type;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			using ancestor::divide_by;
			using complex_toolbox::divide_by;
			NUMERICAL_CONTAINER_CTORS(complex);
			COMPLEX_NUMERICAL_CONTAINER_CTORS;
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(complex<double>(m_value.real().get_d(), m_value.imag().get_d()));
			}
			template <class ArgsTuple>
			complex<double> eval(const double &, const ArgsTuple &) const {
				return complex<double>(m_value.real().get_d(), m_value.imag().get_d());
			}
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
			complex pow(const double &y, const ArgsTuple &) const {
				complex retval;
				retval.m_value = complex<mpq_class>(0, 0);
				// If value = 1, then any power is ok, just return 1.
				if (m_value.real() == 1 and m_value.imag() == 0) {
					retval.m_value.real() = 1;
				} else {
					const piranha::max_fast_int pow_n((piranha::max_fast_int)nearbyint(y));
					if (std::abs(pow_n - y) > piranha::settings::numerical_zero()) {
						// Let's accept real power only of zero value.
						if (m_value.real() != 0 or m_value.imag() != 0) {
							throw(piranha::unsuitable("Cannot raise complex rational coefficient different from unity to real power."));
						}
					} else {
						size_t count;
						if (pow_n >= 0) {
							count = (size_t)pow_n;
							retval.m_value.real() = 1;
							for (size_t i = 0; i < count; ++i) {
								retval.m_value *= m_value;
							}
						} else {
							retval.m_value = m_value;
							const mpq_class div = m_value.real() * m_value.real() + m_value.imag() * m_value.imag();
							retval.m_value.real() /= div;
							retval.m_value.imag() /= div;
							retval.m_value.imag() *= -1;
							count = (size_t)(-pow_n) - 1;
							complex<mpq_class> tmp;
							if (count != 0) {
								tmp = retval.m_value;
							}
							for (size_t i = 0; i < count; ++i) {
								retval.m_value *= tmp;
							}
						}
					}
				}
				return retval;
			}
	};
}

#endif
