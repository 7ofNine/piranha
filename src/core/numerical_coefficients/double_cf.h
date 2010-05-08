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
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <cmath>
#include <complex>

#include "../base_classes/numerical_container.h"
#include "../base_classes/numerical_container_mp.h"
#include "../base_classes/numerical_container_complex_toolbox.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../psym.h"
#include "../settings.h" // Numerical zero.

namespace piranha
{
	// double and complex double will not interoperate with mp types. These template specialisation
	// will convert mp types to double on the fly.
	template <>
	struct numerical_container_value_division<double,mp_rational>
	{
		static void run(double &x, const mp_rational &q)
		{
			x /= q.to_double();
		}
	};

	template <>
	struct numerical_container_value_division<double,mp_integer>
	{
		static void run(double &x, const mp_integer &z)
		{
			x /= z.to_double();
		}
	};

	template <>
	struct numerical_container_value_multiplication<double,mp_rational>
	{
		static void run(double &x, const mp_rational &q)
		{
			x *= q.to_double();
		}
	};

	template <>
	struct numerical_container_value_multiplication<double,mp_integer>
	{
		static void run(double &x, const mp_integer &z)
		{
			x *= z.to_double();
		}
	};

	template <>
	struct numerical_container_value_division<std::complex<double>,mp_rational>
	{
		static void run(std::complex<double> &c, const mp_rational &q)
		{
			c /= q.to_double();
		}
	};

	template <>
	struct numerical_container_value_division<std::complex<double>,mp_integer>
	{
		static void run(std::complex<double> &c, const mp_integer &z)
		{
			c /= z.to_double();
		}
	};

	template <>
	struct numerical_container_value_multiplication<std::complex<double>,mp_rational>
	{
		static void run(std::complex<double> &c, const mp_rational &q)
		{
			c *= q.to_double();
		}
	};

	template <>
	struct numerical_container_value_multiplication<std::complex<double>,mp_integer>
	{
		static void run(std::complex<double> &c, const mp_integer &z)
		{
			c *= z.to_double();
		}
	};

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
			explicit double_cf(): ancestor() {}
			template <class T, class ArgsTuple>
			explicit double_cf(const T &x, const ArgsTuple &args_tuple): ancestor(x,args_tuple) {}
			template <class ArgsTuple>
			explicit double_cf(const mp_rational &q, const ArgsTuple &args_tuple): ancestor(q.to_double(),args_tuple) {}
			template <class ArgsTuple>
			explicit double_cf(const mp_integer &z, const ArgsTuple &args_tuple): ancestor(z.to_double(),args_tuple) {}
			template <class ArgsTuple>
			explicit double_cf(const psym &p, const int &n, const ArgsTuple &a): ancestor(p,n,a) {}
			template <class ArgsTuple>
			double eval(const double &, const ArgsTuple &) const {
				return get_value();
			}
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const
			{
				return std::abs(get_value());
			}
			int to_int() const {
				if (!is_integer(get_value())) {
					piranha_throw(value_error,"cannot convert double coefficient to integer");
				}
				return (int)get_value();
			}
			template <class ArgsTuple>
			std::complex<double_cf> ei(const ArgsTuple &) const;
	};
}

namespace std
{
	template <>
	class complex<piranha::double_cf>:
	public piranha::numerical_container<complex<double>, complex<piranha::double_cf> >,
		public piranha::numerical_container_complex_toolbox<piranha::double_cf>
	{
			typedef piranha::numerical_container<complex<double>, complex<piranha::double_cf> > ancestor;
			typedef piranha::numerical_container_complex_toolbox<piranha::double_cf> complex_toolbox;
		public:
			explicit complex(): ancestor() {}
			template <class T, class ArgsTuple>
			explicit complex(const T &x, const ArgsTuple &args_tuple): ancestor(x,args_tuple) {}
			template <class ArgsTuple>
			explicit complex(const piranha::mp_rational &q, const ArgsTuple &args_tuple): ancestor(q.to_double(),args_tuple) {}
			template <class ArgsTuple>
			explicit complex(const piranha::mp_integer &z, const ArgsTuple &args_tuple): ancestor(z.to_double(),args_tuple) {}
			template <class ArgsTuple>
			explicit complex(const piranha::psym &p, const int &n, const ArgsTuple &a): ancestor(p,n,a) {}
			template <class ArgsTuple>
			complex<double> eval(const double &, const ArgsTuple &) const {
				return get_value();
			}
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const {
				return abs(get_value());
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
		retval.set_value(std::polar(1.,get_value()));
		return retval;
	}
}

#endif
