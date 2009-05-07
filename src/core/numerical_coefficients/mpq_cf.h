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

#include "../base_classes/numerical_container.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"
#include "../settings.h"
#include "../type_traits.h" // For lightweight attribute.
#include "../utils.h" // For is_integer.

namespace piranha
{
	/// Rational numerical coefficient for series.
	/**
	 * Arbitrary-size rational coefficient type, to be used as coefficient in piranha::base_series.
	 */
	class mpq_cf: public numerical_container<mp_rational, mpq_cf>
	{
			// Alias for the parent class.
			typedef numerical_container<mp_rational, mpq_cf> ancestor;
		public:
			typedef mp_rational numerical_type;
			NUMERICAL_CONTAINER_CTORS(mpq_cf)
			// Override norm and evaluation.
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs((double)m_value);
			}
			template <class ArgsTuple>
			double eval(const double &, const ArgsTuple &) const {
				return (double)m_value;
			}
			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (m_value == 0);
			}
			int get_int() const {
				return (int)m_value;
			}
			template <class ArgsTuple>
			mpq_cf pow(const double &y, const ArgsTuple &) const {
				mpq_cf retval;
				retval.m_value = m_value.pow(y);
				return retval;
			}
			template <class ArgsTuple>
			mpq_cf root(const int &n, const ArgsTuple &) const {
				mpq_cf retval;
				retval.m_value = m_value.root(n);
				return retval;
			}
			template <class ArgsTuple>
			std::complex<mpq_cf> ei(const ArgsTuple &) const;
	};
}

namespace std
{
	template <>
	class complex<piranha::mpq_cf>:
				public piranha::numerical_container<std::complex<piranha::mp_rational>, complex<piranha::mpq_cf> >,
				public piranha::numerical_container_complex_toolbox<piranha::mpq_cf>
	{
			typedef piranha::numerical_container<std::complex<piranha::mp_rational>, complex<piranha::mpq_cf> > ancestor;
			typedef piranha::numerical_container_complex_toolbox<piranha::mpq_cf> complex_toolbox;
			friend class piranha::numerical_container_complex_toolbox<piranha::mpq_cf>;
		public:
			using complex_toolbox::divide_by;
			using ancestor::mult_by;
			using complex_toolbox::mult_by;
			typedef piranha::mpq_cf value_type;
			NUMERICAL_CONTAINER_CTORS(complex)
			COMPLEX_NUMERICAL_CONTAINER_CTORS()
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return std::abs((complex<double>)m_value);
			}
			template <class ArgsTuple>
			complex<double> eval(const double &, const ArgsTuple &) const {
				return (complex<double>)(m_value);
			}
			// Override this, hence avoiding to calculate norm.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (m_value == 0);
			}
			template <class ArgsTuple>
			complex pow(const double &y, const ArgsTuple &) const {
				complex retval;
				retval.m_value = m_value.pow(y);
				return retval;
			}
			template <class ArgsTuple>
			complex root(const int &n, const ArgsTuple &) const {
				complex retval;
				retval.m_value =  m_value.root(n);
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
