/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redis\bute it and/or modify  *
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

#ifndef PIRANHA_COMPLEX_GENERIC_MP_CONTAINER_H
#define PIRANHA_COMPLEX_GENERIC_MP_CONTAINER_H

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class T>
	struct cgmp_helper {
		template <class U>
		static void construct(U &u, const T &x)
		{
			u.m_real = x;
		}
		template <class U>
		static void assign(U &u, const T &x)
		{
			u.m_real = x;
			u.m_imag = 0;
		}
		template <class U>
		static bool compare(const U &u, const T &x)
		{
			return (u.m_real == x && u.m_imag == 0);
		}
		template <class U>
		static void add(U &u, const T &x)
		{
			u.m_real += x;
		}
		template <class U>
		static void subtract(U &u, const T &x)
		{
			u.m_real -= x;
		}
		template <class U>
		static void mult(U &u, const T &x)
		{
			u.m_real *= x;
			u.m_imag *= x;
		}
		template <class U>
		static void divide(U &u, const T &x)
		{
			u.m_real /= x;
			u.m_imag /= x;
		}
	};

	template <class T>
	struct cgmp_helper<std::complex<T> > {
		template <class U>
		static void construct(U &u, const std::complex<T> &c)
		{
			u.m_real = c.real();
			u.m_imag = c.imag();
		}
		template <class U>
		static void assign(U &u, const std::complex<T> &c)
		{
			u.m_real = c.real();
			u.m_imag = c.imag();
		}
		template <class U>
		static bool compare(const U &u, const std::complex<T> &c)
		{
			return (u.m_real == c.real() && u.m_imag == c.imag());
		}
		template <class U>
		static void add(U &u, const std::complex<T> &c)
		{
			u.m_real += c.real();
			u.m_imag += c.imag();
		}
		template <class U>
		static void subtract(U &u, const std::complex<T> &c)
		{
			u.m_real -= c.real();
			u.m_imag -= c.imag();
		}
		template <class U>
		static void mult(U &u, const std::complex<T> &c)
		{
			u.mult_by_complex(c);
		}
		template <class U>
		static void divide(U &u, const std::complex<T> &c)
		{
			u.divide_by_complex(c);
		}
	};

	/// Generic container helper for complex multi-precision classes.
	/**
	 * Can be used to help building the complex counterparts of wrapped GMP classes.
	 * The constructors will forward the construction to the
	 * underlying real class. Uses the CRTP.
	 */
	template <class T, class Derived>
	class complex_generic_mp_container
	{
			template <class U>
			friend struct cgmp_helper;
		public:
			/// Get const reference to the real part.
			const T &real() const
			{
				return m_real;
			}
			/// Get const reference to the imaginary part.
			const T &imag() const
			{
				return m_imag;
			}
			/// Set real part.
			void set_real(const T &real) const
			{
				m_real = real;
			}
			/// Set imaginary part.
			void set_imag(const T &imag) const
			{
				m_imag = imag;
			}
			/// Mathematical inversion.
			Derived invert() const
			{
				Derived retval(*derived_const_cast);
				retval.m_imag.negate();
				const T div = abs2();
				retval.m_real /= div;
				retval.m_imag /= div;
				return retval;
			}
			/// Swap content.
			/**
			 * Will call std::swap on real and imaginary members.
			 */
			void swap(Derived &other)
			{
				std::swap(m_real,other.m_real);
				std::swap(m_imag,other.m_imag);
			}
			/// Hash value.
			/**
			 * Uses boost::hash_combine to combine the hashes of real and imaginary parts.
			 */
			size_t hash() const
			{
				boost::hash<T> hasher;
				size_t retval = hasher(m_real);
				boost::hash_combine(retval, hasher(m_imag));
				return retval;
			}
			/// Cast to complex double.
			/**
			 * Will call the to_double() method on real and imaginary part.
			 */
			std::complex<double> to_complex_double() const
			{
				return std::complex<double>(m_real.to_double(),m_imag.to_double());
			}
			/// Cast to complex int.
			/**
			 * Will call the to_int() method on real and imaginary part.
			 */
			std::complex<int> to_complex_int() const
			{
				return std::complex<int>(m_real.to_int(),m_imag.to_int());
			}
		protected:
			/// Default constructor.
			explicit complex_generic_mp_container(): m_real(0),m_imag(0) {}
			/// Generic constructor from single argument.
			/**
			 * Argument can be either a scalar or a complex.
			 */
			template <class U>
			explicit complex_generic_mp_container(const U &x): m_real(0),m_imag(0)
			{
				cgmp_helper<U>::construct(*this,x);
			}
			/// Generic constructor from real and imaginary parts.
			template <class Scalar>
			explicit complex_generic_mp_container(const Scalar &r, const Scalar &i):
				m_real(r),m_imag(i) {}
			/// Constructor from std::string.
			/**
			 * Will raise a value_error exception if string is not valid. A valid string is of the form
			 * "(real,imag)", where "real" and "imag" are strings that can construct successfully an
			 * instance of type T.
			 * @throws value_error if string is invalid.
			 */
			explicit complex_generic_mp_container(const std::string &str): m_real(),m_imag()
			{
				construct_from_string(str.c_str());
			}
			/// Constructor from C string.
			/**
			 * @see complex_generic_mp_container(const std::string &).
			 */
			explicit complex_generic_mp_container(const char *str): m_real(),m_imag()
			{
				construct_from_string(str);
			}
			template <class U>
			Derived &operator=(const U &other)
			{
				cgmp_helper<U>::assign(*this,other);
				return *derived_cast;
			}
			template <class U>
			bool operator==(const U &other) const
			{
				return cgmp_helper<U>::compare(*this,other);
			}
			template <class U>
			Derived &operator+=(const U &other)
			{
				cgmp_helper<U>::add(*this,other);
				return *derived_cast;
			}
			template <class U>
			Derived &operator-=(const U &other)
			{
				cgmp_helper<U>::subtract(*this,other);
				return *derived_cast;
			}
			template <class U>
			Derived &operator*=(const U &other)
			{
				cgmp_helper<U>::mult(*this,other);
				return *derived_cast;
			}
			template <class U>
			Derived &operator/=(const U &other)
			{
				cgmp_helper<U>::divide(*this,other);
				return *derived_cast;
			}
		private:
			// Square of absolute value.
			T abs2() const
			{
				// NOTE: rewrite this in terms of multadd, when implemented.
				T retval(m_real);
				retval *= m_real;
				retval += m_imag * m_imag;
				return retval;
			}
			template <class Complex>
			Derived &divide_by_complex(const Complex &other)
			{
				// NOTE: rewrite this in terms of multadd, when implemented.
				if (other.real() == 0 && other.imag() == 0) {
					piranha_throw(zero_division_error,"cannot divide by zero");
				}
				// This is the divisor, i.e. the square of absolute value of other.
				T div(other.real());
				div *= other.real();
				div += other.imag() * other.imag();
				// The numerator looks like a multiplication with opposite signs.
				const T tmp1(m_imag * other.imag()), tmp2(m_real * other.imag());
				m_imag *= other.real();
				m_imag -= tmp2;
				m_real *= other.real();
				m_real += tmp1;
				// Now divide by divisor.
				m_real /= div;
				m_imag /= div;
				return *derived_cast;
			}
			template <class Complex>
			Derived &mult_by_complex(const Complex &other)
			{
				// NOTE: rewrite this in terms of multadd, when implemented?
				const T tmp1(m_imag * other.imag()), tmp2(m_real * other.imag());
				// NOTE: we do imag first because if we modify real now, then it screws up the computation.
				// m_imag is not used anymore as rhs from this point onwards.
				m_imag *= other.real();
				m_imag += tmp2;
				m_real *= other.real();
				m_real -= tmp1;
				return *derived_cast;
			}
			void construct_from_string(const char *str)
			{
				std::string tmp(str);
				// First let's trim the input string.
				boost::trim(tmp);
				// Next we split the input string into two parts, separated by the comma.
				std::vector<std::string> split_v;
				boost::split(split_v,tmp,boost::is_any_of(","));
				// Now we check that the input string is well formed, i.e., it contains one comma,
				// it starts and ends with round brackets and each element of the split vector contains at
				// least two characters.
				if (split_v.size() != 2 || split_v[0].size() < 2 || split_v[1].size() < 2 ||
					split_v[0][0] != '(' || split_v[1][split_v[1].size() - 1] != ')') {
					piranha_throw(value_error,"invalid string input");
				}
				// Finally, we try to build the two components from the split vector,
				// discarding the positions in which there are the brackets.
				m_real = T(std::string(&split_v[0][1],&split_v[0][split_v[0].size()]));
				m_imag = T(std::string(&split_v[1][0],&split_v[1][split_v[1].size() - 1]));
			}
		protected:
			T	m_real;
			T	m_imag;
	};

	/// Overload out stream operator<<
	template <class T, class Derived>
	inline std::ostream &operator<<(std::ostream &o, const complex_generic_mp_container<T,Derived> &c)
	{
		o << '(' << c.real() << ',' << c.imag() << ')';
		return o;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
