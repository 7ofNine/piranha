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

#ifndef PIRANHA_WIGNER_ROTATION_TOOLBOX_H
#define PIRANHA_WIGNER_ROTATION_TOOLBOX_H

#include <algorithm>
#include <cmath>
#include <complex>

#include "../int_power_cache.h"
#include "../integer_typedefs.h"
#include "../math.h"
#include "../p_assert.h"

namespace piranha
{
	template <class Derived>
	class wigner_rotation_toolbox
	{
		public:
			static std::complex<Derived> Ynm(const max_fast_int &n_, const max_fast_int &m_, const Derived &theta,
				const Derived &phi, const Derived &alpha, const Derived &beta, const Derived &gamma) {
				// Let's fix negative n and/or m.
				max_fast_int n(n_), m(std::abs(m_));
				std::complex<Derived> retval(std::complex<max_fast_int>(static_cast<max_fast_int>(1),
					static_cast<max_fast_int>(0)));
				if (n_ < 0) {
					n = -n_-1;
				}
				if (n == 0 && m == 0) {
					return retval;
				}
				retval = std::complex<Derived>();
				if (m > n) {
					return retval;
				}
				p_assert(n >= m && n >= 0 && m >= 0);
				// Let's prepare the quantities needed for the calculations.
				const std::complex<Derived>
					eit(theta.ei()),
					eib2((beta/static_cast<max_fast_int>(2)).ei());
				const Derived cos_t(eit.real()), sin_t(eit.imag());
				std::complex<Derived> final_factor((gamma*(-m)).ei());
				typedef int_power_cache<Derived> real_cache_type;
				typedef int_power_cache<std::complex<Derived> > complex_cache_type;
				complex_cache_type
					cp(phi.ei(),(phi*static_cast<max_fast_int>(-1)).ei()),
					ca(alpha.ei(),(alpha*static_cast<max_fast_int>(-1)).ei());
				real_cache_type
					ccb2(eib2.real()),
					csb2(eib2.imag());
				final_factor *= einpi2(-m);
				mult_by_factorial(final_factor,n+m);
				for (max_fast_int k = -n; k <= n; ++k) {
					std::complex<Derived> tmp(ca[-k]);
					tmp *= einpi2(k);
					mult_by_factorial(tmp,n-k);
					tmp *= cp[k];
					tmp *= cos_t.Pnm(n,k,sin_t);
					Derived tmp2;
					for (max_fast_int t = std::max<max_fast_int>(static_cast<max_fast_int>(0),k-m);
						t <= std::min<max_fast_int>(n-m,n+k); ++t) {
						Derived tmp3(ccb2[n*2-m+k-t*2]);
						tmp3 *= csb2[m-k+t*2];
						tmp3 *= cs_phase(t);
						div_by_factorial(tmp3,t);
						div_by_factorial(tmp3,n+k-t);
						div_by_factorial(tmp3,n-m-t);
						div_by_factorial(tmp3,m-k+t);
						tmp2 += tmp3;
					}
					tmp *= tmp2;
					retval += tmp;
				}
				retval *= final_factor;
				return retval;
			}
		private:
			template <class T>
			static void mult_by_factorial(T &x, const max_fast_int &n) {
				p_assert(n >= 0);
				for (max_fast_int i = 2; i <= n; ++i) {
					x *= i;
				}
			}
			template <class T>
			static void div_by_factorial(T &x, const max_fast_int &n) {
				p_assert(n >= 0);
				for (max_fast_int i = 2; i <= n; ++i) {
					x /= i;
				}
			}
	};
}

#endif
