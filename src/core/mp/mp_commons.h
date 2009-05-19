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

#ifndef PIRANHA_MP_COMMONS_H
#define PIRANHA_MP_COMMONS_H

namespace piranha
{
	#define derived_const_cast static_cast<Derived const *>(this)
	#define derived_cast static_cast<Derived *>(this)

	/// Toolbox of useful functions for wrapping GMP-like classes.
	/**
	 * CRTP-based.
	 */
	template <class T, class Derived>
	class gmp_toolbox {
		public:
			/// Absolute value.
			Derived abs() const
			{
				return (((*derived_const_cast) >= 0) ? (*derived_const_cast) : -(*derived_const_cast));
			}
			/// Const reference to internal type.
			const T &get_internal() const
			{
				return derived_const_cast->m_value;
			}
			/// Print to stream.
			std::ostream &print(std::ostream &s) const
			{
				return (s << derived_const_cast->m_value);
			}
		protected:
			template <class U>
			bool complex_comparison(const std::complex<U> &c) const
			{
				return (derived_const_cast->m_value == c.real() && c.imag() == 0);
			}
	};

	#undef derived_const_cast
	#undef derived_cast

	// Function for generic exponentiation to rational.
	template <class T, class Rational>
	inline T rat_pow(const T &x, const Rational &q)
	{
		const int n = q.get_num().to_int();
		return x.pow(n).root(q.get_den().to_int());
	}
}

#endif
