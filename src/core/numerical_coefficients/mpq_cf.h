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

#include <boost/type_traits/integral_constant.hpp> // For lightweight attribute.
#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../settings.h"
#include "../type_traits.h" // For lightweight attribute.
#include "../utils.h" // For is_integer.

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
			typedef mpq_class numerical_type;
			// Ctors. Do not use macro since we need to canonicalize.
			/// Empty constructor.
			explicit mpq_cf(): ancestor() {
				m_value.canonicalize();
			}
			/// Constructor from string.
			template <class ArgsTuple>
			explicit mpq_cf(const std::string &s, const ArgsTuple &a): ancestor::numerical_container(s, a) {
				// We need to canonicalize when reading from string.
				m_value.canonicalize();
			}
			/// Constructor from double.
			template <class ArgsTuple>
			explicit mpq_cf(const double &val, const ArgsTuple &a): ancestor::numerical_container(val, a) {
				m_value.canonicalize();
			}
			/// Constructor from psym.
			template <class ArgsTuple>
			explicit mpq_cf(const psym &p, const int &n, const ArgsTuple &a): ancestor::numerical_container(p, n, a) {
				m_value.canonicalize();
			}
			// Override norm and evaluation.
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(m_value.get_d());
			}
			template <class ArgsTuple>
			double eval(const double &, const ArgsTuple &) const {
				return m_value.get_d();
			}
			// Override swapping.
			void swap(mpq_cf &other) {
				mpz_swap(mpq_numref(m_value.get_mpq_t()),mpq_numref(other.m_value.get_mpq_t()));
				mpz_swap(mpq_denref(m_value.get_mpq_t()),mpq_denref(other.m_value.get_mpq_t()));
			}
			// Override division to catch divide by zero.
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
			int get_int() const {
				int retval = m_value.get_num().get_si();
				if (m_value.get_den() != 1) {
					throw(unsuitable("Cannot convert rational coefficient to integer."));
				}
				return retval;
			}
			template <class ArgsTuple>
			mpq_cf pow(const double &y, const ArgsTuple &) const {
				if (utils::is_integer(y)) {
					return pow_int((int)y);
				} else {
				}	return pow_double(y);
			}
			template <class ArgsTuple>
			mpq_cf root(const int &n_, const ArgsTuple &) const {
				mpq_cf retval;
				if (n_ == 0) {
					throw division_by_zero();
				} else if (n_ == 1) {
					retval = *this;
					return retval;
				}
				const size_t n = (n_ > 0) ? n_ : -n_;
				if (!mpz_root(mpq_numref(retval.m_value.get_mpq_t()),mpq_numref(m_value.get_mpq_t()),n) ||
					!mpz_root(mpq_denref(retval.m_value.get_mpq_t()),mpq_denref(m_value.get_mpq_t()),n)) {
					throw unsuitable("Rational coefficient is not an exact nth root.");
				}
				// Better to canonicalise, for peace of mind.
				retval.m_value.canonicalize();
				if (n_ < 0) {
					// Let's guard against division by zero below.
					if (retval.m_value == 0) {
						throw division_by_zero();
					}
					mpq_inv(retval.m_value.get_mpq_t(), mpq_class(retval.m_value).get_mpq_t());
				}
				return retval;
			}
			template <class ArgsTuple>
			std::complex<mpq_cf> ei(const ArgsTuple &) const;
		private:
			mpq_cf pow_int(const int &n) const {
				mpq_cf retval;
				if (m_value == 0 && n < 0) {
					throw division_by_zero();
				}
				if (n < 0) {
					mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (size_t)(-n));
					mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (size_t)(-n));
					// We need to canonicalize, since negative numbers may have gone to the denominator.
					retval.m_value.canonicalize();
				} else {
					mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (size_t)n);
					mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (size_t)n);
				}
				return retval;
			}
			mpq_cf pow_double(const double &y) const {
				mpq_cf retval;
				// If negative, only 1^-something is reasonable.
				if (y < 0) {
					if (m_value == 0) {
						throw division_by_zero();
					} else if (m_value == 1) {
						retval.m_value = 1;
					} else {
						throw unsuitable("Cannot raise rational coefficient different from unity to "
							"negative real power.");
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
						throw unsuitable("Cannot raise rational coefficient different from unity to "
							"positive real power.");
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
			using complex_toolbox::divide_by;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			using complex_toolbox::operator==;
			using complex_toolbox::print_pretty;
			typedef piranha::mpq_cf value_type;
			/// Empty constructor.
			explicit complex(): ancestor() {
				canonicalize();
			}
			/// Constructor from string.
			template <class ArgsTuple>
			explicit complex(const std::string &s, const ArgsTuple &a): ancestor::numerical_container(s, a) {
				// We need to canonicalize when reading from string.
				canonicalize();
			}
			/// Constructor from double.
			template <class ArgsTuple>
			explicit complex(const double &val, const ArgsTuple &a): ancestor::numerical_container(val, a) {
				canonicalize();
			}
			/// Constructor from psym.
			template <class ArgsTuple>
			explicit complex(const piranha::psym &p, const int &n, const ArgsTuple &a):
				ancestor::numerical_container(p, n, a) {
				canonicalize();
			}
			template <class ArgsTuple>
			explicit complex(const std::complex<double> &c, const ArgsTuple &):
				complex_toolbox::numerical_container_complex_toolbox(c) {
				canonicalize();
			}
			template <class ArgsTuple>
			explicit complex(const value_type &r, const ArgsTuple &):
				complex_toolbox::numerical_container_complex_toolbox(r) {
				canonicalize();
			}
			template <class ArgsTuple>
			explicit complex(const value_type &r, const value_type &i, const ArgsTuple &):
				complex_toolbox::numerical_container_complex_toolbox(r, i) {
				canonicalize();
			}
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs(complex<double>(m_value.real().get_d(), m_value.imag().get_d()));
			}
			template <class ArgsTuple>
			complex<double> eval(const double &, const ArgsTuple &) const {
				return complex<double>(m_value.real().get_d(), m_value.imag().get_d());
			}
			// Override swapping.
			void swap(complex &other) {
				mpz_swap(mpq_numref(m_value.real().get_mpq_t()),mpq_numref(other.m_value.real().get_mpq_t()));
				mpz_swap(mpq_denref(m_value.real().get_mpq_t()),mpq_denref(other.m_value.real().get_mpq_t()));
				mpz_swap(mpq_numref(m_value.imag().get_mpq_t()),mpq_numref(other.m_value.imag().get_mpq_t()));
				mpz_swap(mpq_denref(m_value.imag().get_mpq_t()),mpq_denref(other.m_value.imag().get_mpq_t()));
			}
			// Override division to catch divide by zero.
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
				return (m_value.real() == 0 && m_value.imag() == 0);
			}
			template <class ArgsTuple>
			complex pow(const double &y, const ArgsTuple &) const {
				if (piranha::utils::is_integer(y)) {
					return pow_int((int)y);
				} else {
				}	return pow_double(y);
			}
		private:
			void canonicalize() {
				m_value.real().canonicalize();
				m_value.imag().canonicalize();
			}
			complex pow_int(const int &n) const {
				complex retval;
				// For negative powers, we must guard against division by zero.
				if (n < 0) {
					if (m_value.real() == 0 && m_value.imag() == 0) {
						throw piranha::division_by_zero();
						// If source is non-zero, we can invert it and the calculate the power simply by multiplying.
					} else {
						// Invert the source.
						retval.m_value = m_value;
						// div is the module of m_value.
						mpq_class div(m_value.real());
						div *= m_value.real();
						{
							mpq_class tmp(m_value.imag());
							tmp *= m_value.imag();
							div += tmp;
						}
						retval.m_value /= div;
						retval.m_value.imag() *= -1;
						const size_t count = static_cast<size_t>(-n) - 1;
						complex<mpq_class> tmp;
						if (count != 0) {
							tmp = retval.m_value;
						}
						for (size_t i = 0; i < count; ++i) {
							retval.m_value *= tmp;
						}
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
			complex pow_double(const double &y) const {
				complex retval;
				if (y < 0) {
					if (m_value.real() == 0 && m_value.imag() == 0) {
						throw piranha::division_by_zero();
					} else if (m_value.real() == 1 && m_value.imag() == 0) {
						retval.m_value.real() = 1;
						retval.m_value.imag() = 0;
					} else {
						throw piranha::unsuitable("Cannot raise complex rational coefficient "
							"different from unity to negative real power.");
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
						throw piranha::unsuitable("Cannot raise complex rational coefficient "
							"different from unity to positive real power.");
					}
				}
				return retval;
			}
	};
}

namespace piranha
{
	template <class ArgsTuple>
	inline std::complex<mpq_cf> mpq_cf::ei(const ArgsTuple &) const
	{
		throw unsuitable("Rational coefficient is unsuitable for complex exponentiation.");
	}

	// Specialise lightweight type trait for mpq coefficients.
	template <>
	struct is_lightweight<mpq_cf>: boost::true_type {};

	template <>
	struct is_lightweight<std::complex<mpq_cf> >: boost::true_type {};
}

#endif
