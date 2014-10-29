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

#ifndef PIRANHA_COMMON_FOURIER_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_FOURIER_SERIES_TOOLBOX_H

#include <algorithm> // For sorting.
#include <complex>
#include <string>
#include <vector>

#include "../utils.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class common_fourier_series
	{
		public:
			std::complex<Derived> ei() const
			{
				std::complex<Derived> retval(base_ei(derived_const_cast->arguments()));
				retval.setArguments(derived_const_cast->arguments());
				retval.trim();
				return retval;
			}
			Derived cos() const
			{
				return ei().real();
			}
			Derived sin() const
			{
				return ei().imag();
			}
		//protected:
			template <class ArgsTuple>
			std::complex<Derived> base_ei(const ArgsTuple &argsTuple) const
			{
				typedef typename std::complex<Derived>::TermType complex_term_type;
				typedef typename complex_term_type::KeyType KeyType;
				typedef typename Derived::TermType term_type;
				typedef typename std::vector<term_type const *>::const_iterator const_iterator;
				std::complex<Derived> retval;

				if (derived_const_cast->isSingleCf())
                {
					retval.insert(complex_term_type(derived_const_cast->begin()->cf.ei(argsTuple), KeyType()), argsTuple);
				} else
                {
					// Cache and sort the terms according to the criterion defined in the truncator.
					std::vector<term_type const *> cache(derived_const_cast->template get_sorted_series<Derived>(argsTuple));
					// Reverse the series, we want to start multiplication from the least significant terms.
					std::reverse(cache.begin(),cache.end());
					// Let's find out if there is a constant term. If there is one, it will be skipped
					// and multiplied by the result of the Jacobi-Anger expansion of the other terms later.
					// We treat it this way because the constant term may be a phase with arbitrary value,
					// whereas jacang works for coefficients < 1 or so. Since a constant term has all trig
					// multipliers equal to zero, we can just compute its complex exponential without worrying
					// about convergence.
					const_iterator it = cache.begin();
					const const_iterator it_f = cache.end();
					for (; it != it_f; ++it)
                    {
						if ((*it)->key.isUnity())
                        {
							break;
						}
					}
					// Expand using Jacobi-Anger's identity.
					derived_const_cast->jacang(cache, it, retval, argsTuple);
					if (it != it_f)
                    {
						// Take care of the constant element.
						std::complex<Derived> tmp;
						tmp.insert(complex_term_type((*it)->cf.ei(argsTuple),KeyType()),argsTuple);
						retval.baseMultBy(tmp,argsTuple);
					}
				}

				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
