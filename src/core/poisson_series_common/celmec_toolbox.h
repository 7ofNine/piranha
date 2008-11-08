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

#ifndef PIRANHA_CELMEC_TOOLBOX_H
#define PIRANHA_CELMEC_TOOLBOX_H

#include <string>

#include "../integer_typedefs.h"
#include "../psym.h"

namespace piranha
{
	/// Toolbox for Celestial Mechanics.
	/**
	 * This toolbox requires the following toolboxes: multiplication toolbox, Poisson series toolbox,
	 * special functions toolbox. Derived must be a Poisson series with
	 * polynomial arguments in slot 0 of the arguments tuple.
	 */
	template <class Derived>
	class celmec_toolbox
	{
		public:
			static Derived r_a(const Derived &e_series, const Derived &M_series) {
				// First let's build 1+1/2 e^2.
				Derived retval(e_series);
				retval *= e_series;
				retval /= (max_fast_int)2;
				retval += (max_fast_int)1;
				// Let's find out the upper limit of the r_a development, according to the truncation limits
				// set in the truncator.
				const size_t n = e_series.psi();
				Derived tmp;
				for (size_t i = 1; i <= n; ++i) {
					Derived expansion_term(e_series);
					expansion_term *= (max_fast_int)i;
					expansion_term = expansion_term.dbesselJ((max_fast_int)i);
					expansion_term /= (max_fast_int)i;
					Derived trig(M_series);
					trig *= (max_fast_int)i;
					trig = trig.cos();
					expansion_term *= trig;
					tmp += expansion_term;
				}
				tmp *= e_series;
				tmp *= (max_fast_int)(-2);
				retval += tmp;
				return retval;
			}
			static Derived cos_E(const Derived &e_series, const Derived &M_series) {
				// First let's build 1/2 e.
				Derived retval(e_series);
				retval /= (max_fast_int)2;
				const size_t n = e_series.psi();
				Derived tmp;
				for (size_t i = 1; i <= n; ++i) {
					Derived expansion_term(e_series);
					expansion_term *= (max_fast_int)i;
					expansion_term = expansion_term.dbesselJ((max_fast_int)i);
					expansion_term /= (max_fast_int)i;
					Derived trig(M_series);
					trig *= (max_fast_int)i;
					trig = trig.cos();
					expansion_term *= trig;
					tmp += expansion_term;
				}
				tmp *= (max_fast_int)(-2);
				retval += tmp;
				retval *= (max_fast_int)(-1);
				return retval;
			}
			static Derived sin_E(const Derived &e_series, const Derived &M_series) {
				const size_t n = e_series.psi();
				Derived retval;
				for (size_t i = 1; i <= n; ++i) {
					Derived expansion_term(e_series);
					expansion_term *= (max_fast_int)i;
					expansion_term = expansion_term.besselJ_div_m((max_fast_int)i,1);
					Derived trig(M_series);
					trig *= (max_fast_int)i;
					trig = trig.sin();
					expansion_term *= trig;
					retval += expansion_term;
				}
				retval *= (max_fast_int)(2);
				return retval;
			}
			static Derived sin_f(const Derived &e_series, const Derived &M_series) {
				Derived tmp(e_series);
				tmp *= e_series;
				tmp = (max_fast_int)1 - tmp;
				tmp = tmp.root(2);
				tmp *= (max_fast_int)2;
				Derived retval;
				const size_t n = e_series.psi();
				// Regarding range here: we must perform n iterations for the power series, so starting
				// from i = 1 we must reach n included.
				for (size_t i = 1; i <= n; ++i) {
					Derived expansion_term(e_series);
					expansion_term *= (max_fast_int)i;
					expansion_term = expansion_term.dbesselJ((max_fast_int)i);
					Derived trig(M_series);
					trig *= (max_fast_int)i;
					trig = trig.sin();
					expansion_term *= trig;
					retval += expansion_term;
				}
				retval *= tmp;
				return retval;
			}
			static Derived cos_f(const Derived &e_series, const Derived &M_series) {
				// 2*(1-e**2).
				Derived tmp(e_series);
				tmp *= e_series;
				tmp = (max_fast_int)1 - tmp;
				tmp *= (max_fast_int)2;
				Derived retval;
				const size_t n = e_series.psi();
				// See above.
				for (size_t i = 1; i <= n; ++i) {
					Derived expansion_term(e_series);
					expansion_term *= (max_fast_int)i;
					expansion_term = expansion_term.besselJ_div_m((max_fast_int)i,1);
					expansion_term *= (max_fast_int)i;
					Derived trig(M_series);
					trig *= (max_fast_int)i;
					trig = trig.cos();
					expansion_term *= trig;
					retval += expansion_term;
				}
				retval *= tmp;
				retval -= e_series;
				return retval;
			}
			static Derived EE(const Derived &e_series, const Derived &M_series) {
				Derived retval(M_series);
				const size_t n = e_series.psi(1);
				Derived tmp;
				for (size_t i = 1; i <= n; ++i) {
					Derived expansion_term((e_series * (max_fast_int)i).besselJ((max_fast_int)i));
					expansion_term /= (max_fast_int)i;
					expansion_term *= (M_series * (max_fast_int)i).sin();
					tmp += expansion_term;
				}
				tmp *= (max_fast_int)2;
				retval += tmp;
				return retval;
			}
			static std::complex<Derived> eipE(const Derived &e, const Derived &M, const max_fast_int &p) {
				// S1
				std::complex<Derived> S1;
				if (p > 0) {
					for (max_fast_int n = 1; n < p; ++n) {
						std::complex<Derived> tmp = (M * n).ei();
						tmp *= (e * n).besselJ(p - n);
						tmp /= n;
						tmp *= cs_phase(p - n);
						S1 += tmp;
					}
					const max_fast_int iter = e.psi();
					for (max_fast_int n = p; n < p + iter; ++n) {
						std::complex<Derived> tmp = (M * n).ei();
						tmp *= (e * n).besselJ(n - p);
						tmp /= n;
						S1 += tmp;
					}
				} else {
					const max_fast_int iter = e.psi(1 - p);
					for (max_fast_int n = 1; n <= iter; ++n) {
						std::complex<Derived> tmp = (M * n).ei();
						tmp *= (e * n).besselJ(n - p);
						tmp /= n;
						S1 += tmp;
					}
				}
				S1 *= p;
				// S2
				std::complex<Derived> S2;
				if (p > 0) {
					const max_fast_int iter = e.psi(1 + p);
					for (max_fast_int n = 1; n <= iter; ++n) {
						std::complex<Derived> tmp = (M * -n).ei();
						tmp *= (e * -n).besselJ(n + p);
						tmp *= cs_phase(n + p);
						tmp /= n;
						S2 += tmp;
					}
				} else {
					for (max_fast_int n = 1; n < -p; ++n) {
						std::complex<Derived> tmp = (M * -n).ei();
						tmp *= (e * -n).besselJ(-n - p);
						tmp /= n;
						S2 += tmp;
					}
					const max_fast_int iter = e.psi();
					for (max_fast_int n = -p; n < -p + iter; ++n) {
						std::complex<Derived> tmp = (M * -n).ei();
						tmp *= (e * -n).besselJ(n + p);
						tmp *= cs_phase(n + p);
						tmp /= n;
						S2 += tmp;
					}
				}
				S2 *= -p;
				// c0
				Derived c0;
				switch (p) {
					case 0:
						c0 = Derived(max_fast_int(1));
						break;
					case 1:
						c0 = Derived(e / max_fast_int(-2));
						break;
					case -1:
						c0 = Derived(e / max_fast_int(-2));
						break;
					default:
						;
				}
				// Compose the various parts and return.
				S1 += S2;
				S1 += c0;
				return S1;
			}
	};
}

#endif
