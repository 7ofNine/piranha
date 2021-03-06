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

#ifndef PIRANHA_FOURIER_SERIES_TERM_H
#define PIRANHA_FOURIER_SERIES_TERM_H

#include <boost/tuple/tuple.hpp>

#include "../base_classes/base_term.h"

namespace piranha
{
	// Term class for Fourier series.
	// Cf:   coefficients for series term e.g. double_cf, PolynomialCf
	// Trig: trigonometric key i.e. the linear constituents and their coefficients and cos/sin (flavour), e.g. TrigVector<int16_t, 0>, TrigVectorint16_t, 1>,
	//       last template parameter is actually the echelon level.
	// Separator: print/read separator between coefficient and key, e.g.:  '|'
	// Allocator: specific allocator e.g. for statistics or performance improvements. but typicall std::allocator<char>
    // 
    // there is no verification that Trig is actually a TrigVector class
	template <class Cf, class Trig, char Separator>
	class FourierSeriesTerm: public BaseTerm<Cf, Trig, Separator, FourierSeriesTerm<Cf, Trig, Separator> >
	{
			// Alias for the ancestor.
			typedef BaseTerm<Cf, Trig, Separator, FourierSeriesTerm> ancestor;

		public:

			
			typedef Cf CfType;   // Alias for coefficient type.
			
			typedef Trig KeyType; // Alias for trigonometric type.
			
			typedef typename boost::tuple<FourierSeriesTerm, FourierSeriesTerm> multiplication_result; // Result of the multiplication of two terms. they result is two cos or sin terms

			PIRANHA_TERM_CTORS(FourierSeriesTerm);

			template <class ArgsTuple>                   
			bool isCanonical(const ArgsTuple &) const    // Check if the term is canonical.
            {
				return (ancestor::key.sign() > 0);
			}


			// TODO: check if it makes sense to skip the check here and assume canonicalise will be used iff
			// isCanonical has already been tested.
			
			template <class ArgsTuple>                     
			void canonicalise(ArgsTuple const &argsTuple)    // Canonicalise the term.
            {
				if (!isCanonical(argsTuple)) 
                {
					invertTrigArgs(argsTuple);
				}
			}


			
			//
			// NOTE: the result of multiplication here _must_ be canonical.
			//
			template <class Term1, class Term2, class ArgsTuple>
			static void multiply(const Term1 &t1, const Term2 &t2, multiplication_result & res, const ArgsTuple &argsTuple)  // Term multiplication.
            {
				// trigonometric multiplication.
				t1.key.multiply(t2.key, res.template get<0>().key, res.template get<1>().key);
				
				// coefficient multiplication. First coefficient, then assign the second one.
				// 
				// TODO: maybe provide the semantics to coefficients for something like this:
				// cf1.multiply_by_cf(cf2,res.template get<0>().cf,argsTuple),
				// so that we can avoid a copy.
				res.template get<0>().cf = t1.cf;
				res.template get<0>().cf.multBy(t2.cf, argsTuple);
				res.template get<0>().cf.divideBy(2, argsTuple);
				res.template get<1>().cf = res.template get<0>().cf;

				// adjust the signs according to werner's formulas. 
				if (t1.key.getFlavour() == t2.key.getFlavour()) 
                {
					res.template get<0>().key.setFlavour(true);
					res.template get<1>().key.setFlavour(true);
					if (!t1.key.getFlavour()) 
                    {
						res.template get<1>().cf.invertSign(argsTuple);
					}
				} else 
                {
					res.template get<0>().key.setFlavour(false);
					res.template get<1>().key.setFlavour(false);
					if (t1.key.getFlavour()) 
                    {
						res.template get<0>().cf.invertSign(argsTuple);
					}
				}

				// canonicalise the retval terms.
				res.template get<0>().canonicalise(argsTuple);
				res.template get<1>().canonicalise(argsTuple);
			}

		private:

			// invert sign of term.
			template <class ArgsTuple>
			void invertTrigArgs(const ArgsTuple &argsTuple) 
            {
				ancestor::key.invertSign();
				if (!(ancestor::key.getFlavour())) 
                {
					ancestor::cf.invertSign(argsTuple);
				}
			}
	};
}
#endif
