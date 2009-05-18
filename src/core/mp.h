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

#ifndef PIRANHA_MP_H
#define PIRANHA_MP_H

#include <cmath>
#include <complex>

// For now we support only GMP.
#include "mp/piranha_gmp.h"

namespace piranha
{

}

namespace std
{
	/// Overload std::pow for double and piranha::mp_rational arguments.
	inline double pow(const double &x, const piranha::mp_rational &q)
	{
		return pow(x,q.to_double());
	}

	/// Overload std::pow for std::complex<double> and piranha::mp_rational arguments.
	inline complex<double> pow(const complex<double> &c, const piranha::mp_rational &q)
	{
		return pow(c,q.to_double());
	}

	/// Overload std::pow for int and piranha::mp_rational arguments.
	inline double pow(const int &n, const piranha::mp_rational &q)
	{
		return pow(n,q.to_double());
	}

	/// Overload std::pow for std::complex<int> and piranha::mp_rational arguments.
	inline complex<int> pow(const complex<int> &c, const piranha::mp_rational &q)
	{
		return pow(c,q.to_double());
	}
}

#endif
