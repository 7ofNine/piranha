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

#ifndef PIRANHA_BASE_SERIES_COMPLEX_TOOLBOX_H
#define PIRANHA_BASE_SERIES_COMPLEX_TOOLBOX_H

#include <complex>

#include "../exceptions.h"
#include "base_series.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class BaseSeriesComplex
	{
		public:
			typedef std::complex<RealDerived> Derived;
			typedef RealDerived value_type;

			template <class ArgsTuple>
			RealDerived baseReal(const ArgsTuple &argsTuple) const 
			{
				return getComponent<0>(argsTuple);
			}


			template <class ArgsTuple>
			RealDerived baseImag(const ArgsTuple &argsTuple) const 
			{
				return getComponent<1>(argsTuple);
			}


			template <class ArgsTuple>
			void baseSetReal(const RealDerived &real, const ArgsTuple &argsTuple) 
			{
				derived_cast->baseSubtract(baseReal(argsTuple), argsTuple);
				derived_cast->baseAdd(real, argsTuple);
			}


			template <class ArgsTuple>
			void baseSetImag(const RealDerived &imaginary, const ArgsTuple &argsTuple) 
			{
				typedef typename RealDerived::const_iterator RealIterator;
				typedef typename Derived::TermType ComplexTermType;

				ComplexTermType tmp;
				// First let's remove the old imaginary part.
				RealDerived oldImaginary(baseImag(argsTuple));
				const RealIterator oldImaginaryEnd = oldImaginary.end();

				for (RealIterator it = oldImaginary.begin(); it != oldImaginaryEnd; ++it) 
				{
					tmp.key = it->key;
					tmp.cf.setImag(it->cf, argsTuple);
					derived_cast->template insert<true, false>(tmp, argsTuple);
				}

				// Now add the new imaginary part.
				const RealIterator itEnd = imaginary.end();
				
                for (RealIterator it = imaginary.begin(); it != itEnd; ++it) 
				{
					tmp.key = it->key;
					tmp.cf.setImag(it->cf, argsTuple);
					derived_cast->insert(tmp, argsTuple);
				}
			}


            // square of complex norm/absolute value
			template <class ArgsTuple>
			RealDerived baseAbs2(const ArgsTuple &argsTuple) const
			{
				RealDerived retval = baseReal(argsTuple);
                RealDerived tmp    = baseImag(argsTuple);
				retval.baseMultBy(retval, argsTuple);
				tmp.baseMultBy(tmp, argsTuple);
				retval.baseAdd(tmp, argsTuple);
				return retval;
			}

            // complex norm/absolute value
            template <class ArgsTuple>
			RealDerived baseAbs(const ArgsTuple &argsTuple) const
			{
				return baseAbs2(argsTuple).baseRoot(2, argsTuple);
			}

            // complex conjugate
			template <class ArgsTuple>
			Derived baseConjugate(const ArgsTuple &argsTuple) const 
			{
				Derived retval;
				retval.baseAdd(baseReal(argsTuple), argsTuple);
				RealDerived tmp = baseImag(argsTuple);
				tmp.baseMultBy(-1, argsTuple);
				retval.baseSetImag(tmp, argsTuple);
				return retval;
			}


			// Specialise inversion to use conjugate * inverse of absolute value ** 2. This is useful
			// when the complex series is a complex exponential of something.
			template <class ArgsTuple>
			Derived baseInvert(const ArgsTuple &argsTuple) const 
			{
				Derived retval = baseConjugate(argsTuple);
				retval.baseMultBy(baseAbs2(argsTuple).basePow(-1, argsTuple), argsTuple);
				return retval;
			}


			template <class ArgsTuple>
			void baseConstructFromReal(const RealDerived &real, const ArgsTuple &argsTuple) 
			{
				// Make sure we are being called from an empty series.
				PIRANHA_ASSERT(derived_const_cast->empty());
				derived_cast->insertRange(real.begin(), real.end(), argsTuple);
			}


			template <class ArgsTuple>
			void baseConstructFromRealImag(const RealDerived &real, const RealDerived &imaginary, const ArgsTuple &argsTuple) 
			{
				typedef typename RealDerived::const_iterator RealIterator;
				typedef typename Derived::TermType ComplexTermType;

				// Make sure we are being called from an empty series.
				PIRANHA_ASSERT(derived_const_cast->empty());
				
                // Let's build the real part first.
				baseConstructFromReal(real, argsTuple);
				// Now let's proceed to the imaginary part.
				const RealIterator itEnd = imaginary.end();
				ComplexTermType tmp;

				for (RealIterator it = imaginary.begin(); it != itEnd; ++it) 
				{
					tmp.key = it->key;
					tmp.cf.setImag(it->cf, argsTuple);
					derived_cast->insert(tmp, argsTuple);
				}
			}

		private:
			// Use N = 0 for real, N != 0 for imag.
			template <int N, class Real, class ArgsTuple>
			static Real getCfComponent(const std::complex<Real> &c, const ArgsTuple &argsTuple) 
			{
				if (N) 
				{
					return c.imag(argsTuple);
				} else 
				{
					return c.real(argsTuple);
				}
			}


			template <int N, class ArgsTuple>
			RealDerived getComponent(const ArgsTuple &argsTuple) const 
			{
				typedef typename Derived::const_iterator ComplexIterator;
				RealDerived retval;
				
                const ComplexIterator itEnd = derived_const_cast->end();
				for (ComplexIterator it = derived_const_cast->begin(); it != itEnd; ++it) 
				{
					typename RealDerived::TermType tmp(getCfComponent<N>(it->cf, argsTuple), it->key);
					retval.insert(tmp, argsTuple);
				}

				return retval;
			}
	};

#define COMPLEX_E0_SERIES_TERM(TermName) TermName<std::complex<Cf>, Key, '|', Allocator>
#define COMPLEX_E0_SERIES(SeriesName) std::complex<E0_SERIES(SeriesName)>
#define COMPLEX_E0_SERIES_BASE_ANCESTOR(TermName, SeriesName) piranha::BaseSeries<COMPLEX_E0_SERIES_TERM(TermName), '\n', \
	Allocator, COMPLEX_E0_SERIES(SeriesName) >

#define COMPLEX_E1_SERIES_TERM(TermName, CfName) TermName<std::complex<CfName>, Key1, '|', Allocator>
#define COMPLEX_E1_SERIES(SeriesName) std::complex<E1_SERIES(SeriesName)>
#define COMPLEX_E1_SERIES_BASE_ANCESTOR(TermName, CfName, SeriesName) piranha::BaseSeries<COMPLEX_E1_SERIES_TERM(TermName, CfName), \
	'\n', Allocator,COMPLEX_E1_SERIES(SeriesName) >
}

#undef derived_const_cast
#undef derived_cast

#endif
