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
#include <vector>

#include "../common_functors.h"
#include "../exceptions.h"
#include "../math.h"
#include "../power_cache.h"
#include "../mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class NamedSeriesSpecialFunctions
	{
		public:

			/// Bessel function of the first kind.
			Derived besselJ(const int &order) const 
            {
				Derived retval(derived_const_cast->baseBesselJ(order, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			Derived dbesselJ(const int &order) const 
            {
				Derived retval(derived_const_cast->baseDBesselJ(order, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			/// Bessel function of the first kind of integer order divided by its argument**m.
			Derived besselJ_div_m(const int &order, const int &m) const 
            {
				Derived retval(derived_const_cast->baseBesselJDivm(order, m, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			Derived hyperF(const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list, const int &iter_limit) const
			{
				if (iter_limit < 0)
                {
					PIRANHA_THROW(value_error,"iteration limit in hyperF must be non-negative");
				}

				Derived retval(derived_const_cast->baseHyperF(a_list, b_list, iter_limit, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			Derived hyperF(const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list) const
			{
				Derived retval(derived_const_cast->baseHyperF(a_list, b_list, -1, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			Derived dhyperF(const int &n, const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list, const int &iter_limit) const
			{
				if (iter_limit < 0) {
					PIRANHA_THROW(value_error,"iteration limit in dhyperF must be non-negative");
				}
				Derived retval(derived_const_cast->baseDHyperF(n, a_list, b_list, iter_limit, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			Derived dhyperF(const int &n, const std::vector<mp_rational> &a_list, const std::vector<mp_rational> &b_list) const
			{
				Derived retval(derived_const_cast->baseDHyperF(n, a_list, b_list, -1, derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}


			/// Associated Legendre Function of degree n and order m.
			/**
			 * Implemented through multiple differentiation of hypergeometric series.
			 * When m is odd, this method will attempt to calculate sqrt(1 - self ** 2) internally.
			 */
			Derived legendrePnm(const int &n, const int &m) const
			{
				return impl_legendrePnm(n, m, 0);
			}


			/// Associated Legendre Function of degree n and order m.
			/**
			 * The additional input parameter self_qc is the quadratic conjugate of self, i.e. sqrt(1 - self ** 2).
			 * It will be used when m is odd.
			 */
			Derived legendrePnm(const int &n, const int &m, const Derived &self_qc) const
			{
				return impl_legendrePnm(n, m, &self_qc);
			}


			/// Legendre polynomial.
			/**
			 * Equivalent to legendrePnm(n,0).
			 */
			Derived legendrePn(const int &n) const
			{
				if (n < 0) {
					PIRANHA_THROW(value_error,"please select a non-negative order");
				}
				return legendrePnm(n, 0);
			}


			static std::complex<Derived> Ynm(const int &n, const int &m, const Derived &theta, const Derived &phi)
            {
				const std::complex<Derived> ei_theta(theta.ei());
				std::complex<Derived> retval((phi * m).ei());
				retval *= ei_theta.real().legendrePnm(n, m, ei_theta.imag());
				return retval;
			}


			static std::complex<Derived> Ynm(const int &n, const int &m, const Derived &theta, const std::complex<Derived> &ei_phi, const std::complex<Derived> &emi_phi) 
            {
				const std::complex<Derived> ei_theta(theta.ei());
				std::complex<Derived> retval;

				if (m >= 0) 
                {
					retval = ei_phi.pow(m);
				} else 
                {
					retval = emi_phi.pow(m);
				}

				retval *= ei_theta.real().legendrePnm(n, m, ei_theta.imag());
				return retval;
			}


			static std::complex<Derived> Ynm(const int &n_, const int &m_, const Derived &theta, const std::complex<Derived> &ei_phi, const std::complex<Derived> &emi_phi,
				                             const Derived &alpha, const Derived &beta, const Derived &gamma) 
            {
				// Let's fix negative n and/or m.
				int n(n_), m(std::abs(m_));
				std::complex<Derived> retval(std::complex<double>(1,0));

				if (n_ < 0) 
                {
					n = -n_-1;
				}

				if (n == 0 && m == 0) 
                {
					return retval;
				}

				retval = std::complex<Derived>();
				if (m > n) 
                {
					return retval;
				}

				PIRANHA_ASSERT(n >= m && n >= 0 && m >= 0);

				// Let's prepare the quantities needed for the calculations.
				const std::complex<Derived>
					eit(theta.ei()),
					eib2((beta/2).ei());
				const Derived cos_t(eit.real()), sin_t(eit.imag());
				std::complex<Derived> final_factor((gamma*(-m)).ei());
				typedef PowerCache<Derived, int, NamedSeriesArithmetics<Derived> > real_cache_type;
				typedef PowerCache<std::complex<Derived>,int,
					NamedSeriesArithmetics<std::complex<Derived> > > complex_cache_type;

				complex_cache_type
					cp(ei_phi,emi_phi),
					ca(alpha.ei(),(alpha * -1).ei());
				real_cache_type
					ccb2(eib2.real()),
					csb2(eib2.imag());
				final_factor *= einpi2(-m);
				final_factor *= factorial(n+m);

				for (int k = -n; k <= n; ++k) 
                {
					std::complex<Derived> tmp(ca[-k]);
					tmp *= einpi2(k);
					tmp *= factorial(n-k);
					tmp *= cp[k];
					tmp *= cos_t.legendrePnm(n,k,sin_t);
					Derived tmp2;
					for (int t = std::max<int>(0,k-m); t <= std::min<int>(n-m,n+k); ++t) 
                    {
						Derived tmp3(ccb2[n*2-m+k-t*2]);
						tmp3 *= csb2[m-k+t*2];
						tmp3 *= cs_phase(t);
						tmp3 /= factorial(t);
						tmp3 /= factorial(n+k-t);
						tmp3 /= factorial(n-m-t);
						tmp3 /= factorial(m-k+t);
						tmp2 += tmp3;
					}

					tmp *= tmp2;
					retval += tmp;
				}

				retval *= final_factor;
				return retval;
			}


			static std::complex<Derived> Ynm(const int &n, const int &m, const Derived &theta,
				const Derived &phi, const Derived &alpha, const Derived &beta, const Derived &gamma) 
            {
				return Ynm(n,m,theta,phi.ei(),(phi * -1).ei(),alpha,beta,gamma);
			}


		private:

			Derived impl_legendrePnm(const int &n_, const int &m_, const Derived *self_qc) const
			{
				// Take care of negative degree and/or order.
				const int n = (n_ >= 0) ? n_ : (-n_ - 1), m = (m_ >= 0) ? m_ : -m_;
				Derived retval;
				if (m > n) 
                {
					return retval;
				}

				std::vector<mp_rational> a_list, b_list;
				a_list.push_back(mp_rational(-n));
				a_list.push_back(mp_rational(n) + 1);
				b_list.push_back(mp_rational(1));
				retval = ((1 - *derived_const_cast) / 2).dhyperF(m,a_list,b_list);
				// If m is odd we need to deal with a square root.
				if (m & 1) 
                {
					if (self_qc) 
                    {
						retval *= self_qc->pow(m);
					} else 
                    {
						retval *= (1 - derived_const_cast->pow(2)).root(2).pow(m);
					}
				} else 
                {
					retval *= (1 - derived_const_cast->pow(2)).pow(m / 2);
				}
				// Correct for the fact that we are not deriving with respect to self but to
				// (1 - self) / 2. The cs_phase of the original definition of Pnm and the one
				// of this correction cancel each other out.
				retval /= mp_integer(2).pow(m);
				// Finally, if original m was negative, we must multiply by another correcting factor.
				if (m_ < 0) 
                {
					retval *= cs_phase(m);
					retval *= (mp_integer(n) - m).factorial();
					retval /= (mp_integer(n) + m).factorial();
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
