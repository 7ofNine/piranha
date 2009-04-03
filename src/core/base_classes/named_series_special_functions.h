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

#ifndef PIRANHA_NAMED_SERIES_SPECIAL_FUNCTIONS_H
#define PIRANHA_NAMED_SERIES_SPECIAL_FUNCTIONS_H

#include <algorithm>
#include <cmath>
#include <complex>

#include "../common_functors.h"
#include "../int_power_cache.h"
#include "../math.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class named_series_special_functions
	{
		public:
			/// Bessel function of the first kind.
			Derived besselJ(const int &order) const {
				Derived retval(derived_const_cast->besselJ_(order, derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			Derived dbesselJ(const int &order) const {
				Derived retval(derived_const_cast->dbesselJ_(order, derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			/// Bessel function of the first kind of integer order divided by its argument**m.
			Derived besselJ_div_m(const int &order, const int &m) const {
				Derived retval(derived_const_cast->besselJ_div_m_(order, m, derived_const_cast->m_arguments));
				retval.m_arguments = derived_const_cast->m_arguments;
				retval.trim();
				return retval;
			}
			/// Legendre function of the first kind: Pnm(self).
			/**
			* This implementation uses recurrence relations. self_qc is the quadratic conjugate of self, i.e.,
			* sqrt(1-self**2).
			*/
			// This recursion appears on Wikipedia.
			Derived Pnm(const int &n_, const int &m_, const Derived &self_qc) const {
				// This is P00 right now.
				Derived retval(1);
				int n(n_), m(std::abs(m_));
				if (n_ < 0) {
					n = -n_-1;
				}
				if (n == 0 && m == 0) {
					return retval;
				}
				if (m > n) {
					retval = Derived();
					return retval;
				}
				p_assert(n >= m && n >= 0 && m >= 0);
				Derived P00(retval), old_Pnm, tmp1;
				int i = 0;
				// Recursion to get from P_00 to P_mm.
				for (; i < m; ++i) {
					retval *= -(i*2+1);
					retval *= self_qc;
				}
				p_assert(i == m);
				// Recursion to get from P_mm to P_nm (n>m).
				for (; i < n; ++i) {
					old_Pnm *= -m-i;
					old_Pnm /= -m+i+1;
					tmp1 = old_Pnm;
					old_Pnm = retval;
					retval *= i*2+1;
					retval /= -m+i+1;
					retval *= *derived_const_cast;
					retval += tmp1;
				}
				if (m_ < 0) {
					for (int j = n + m; j >= n - m + 1; --j) {
						retval /= j;
					}
					retval *= cs_phase(m);
				}
				return retval;
			}
			Derived Pnm(const int &n, const int &m) const {
				// Let's build sqrt(1-self**2).
				Derived tmp(*derived_const_cast);
				tmp *= tmp;
				tmp *= -1;
				tmp += Derived(1);
				return Pnm(n,m,tmp.root(2));
			}
			Derived Pn(const int &n) const {
				return Pnm(n,0,Derived());
			}
			static std::complex<Derived> Ynm(const int &n, const int &m,
				const Derived &theta, const Derived &phi) {
				const std::complex<Derived> ei_theta(theta.ei());
				std::complex<Derived> retval((phi * m).ei());
				retval *= ei_theta.real().Pnm(n,m,ei_theta.imag());
				return retval;
			}
			static std::complex<Derived> Ynm(const int &n, const int &m,
				const Derived &theta, const std::complex<Derived> &ei_phi, const std::complex<Derived> &emi_phi) {
				const std::complex<Derived> ei_theta(theta.ei());
				std::complex<Derived> retval;
				if (m >= 0) {
					retval = ei_phi.pow(m);
				} else {
					retval = emi_phi.pow(m);
				}
				retval *= ei_theta.real().Pnm(n,m,ei_theta.imag());
				return retval;
			}
			static std::complex<Derived> Ynm(const int &n_, const int &m_, const Derived &theta,
				const std::complex<Derived> &ei_phi, const std::complex<Derived> &emi_phi,
				const Derived &alpha, const Derived &beta, const Derived &gamma) {
				// Let's fix negative n and/or m.
				int n(n_), m(std::abs(m_));
				std::complex<Derived> retval(std::complex<double>(1,0));
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
					eib2((beta/2).ei());
				const Derived cos_t(eit.real()), sin_t(eit.imag());
				std::complex<Derived> final_factor((gamma*(-m)).ei());
				typedef int_power_cache<Derived,named_series_arithmetics<Derived> > real_cache_type;
				typedef int_power_cache<std::complex<Derived>,
					named_series_arithmetics<std::complex<Derived> > > complex_cache_type;
				complex_cache_type
					cp(ei_phi,emi_phi),
					ca(alpha.ei(),(alpha * -1).ei());
				real_cache_type
					ccb2(eib2.real()),
					csb2(eib2.imag());
				final_factor *= einpi2(-m);
				final_factor *= final_factor.factorial(n+m);
				for (int k = -n; k <= n; ++k) {
					std::complex<Derived> tmp(ca[-k]);
					tmp *= einpi2(k);
					tmp *= tmp.factorial(n-k);
					tmp *= cp[k];
					tmp *= cos_t.Pnm(n,k,sin_t);
					Derived tmp2;
					for (int t = std::max<int>(0,k-m);
						t <= std::min<int>(n-m,n+k); ++t) {
						Derived tmp3(ccb2[n*2-m+k-t*2]);
						tmp3 *= csb2[m-k+t*2];
						tmp3 *= cs_phase(t);
						tmp3 *= tmp3.factorial(t).inv();
						tmp3 *= tmp3.factorial(n+k-t).inv();
						tmp3 *= tmp3.factorial(n-m-t).inv();
						tmp3 *= tmp3.factorial(m-k+t).inv();
						tmp2 += tmp3;
					}
					tmp *= tmp2;
					retval += tmp;
				}
				retval *= final_factor;
				return retval;
			}
			static std::complex<Derived> Ynm(const int &n, const int &m, const Derived &theta,
				const Derived &phi, const Derived &alpha, const Derived &beta, const Derived &gamma) {
				return Ynm(n,m,theta,phi.ei(),(phi * -1).ei(),alpha,beta,gamma);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
