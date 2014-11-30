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
	// Cf:   coefficients for series term e.g. double_cf, polynomial_cf
	// Trig: trigonometric key i.e. te linear constituents and their coefficients and cos/sin (flavour), e.g. TrigVector<boost::int16_t, 0>, TrigVector<boost::int16_t, 1>,
	//       last template parameter is actually the echelon level.
	// Separataor: print/read separator between coefficient and key, e.g.:  '|'
	// Allocator: specific allocator e.g. for statistics or performance improvements. but typicall std::allocator<char>
	template <class Cf, class Trig, char Separator, class Allocator>
	class FourierSeriesTerm: public BaseTerm<Cf, Trig, Separator, Allocator, FourierSeriesTerm<Cf, Trig, Separator, Allocator> >
	{
			// Alias for the ancestor.
			typedef BaseTerm<Cf, Trig, Separator, Allocator, FourierSeriesTerm> ancestor;

		public:

			/// Alias for coefficient type.
			typedef Cf CfType;
			/// Alias for trigonometric type.
			typedef Trig KeyType;
			/// Result of the multiplication of two terms.
			typedef typename boost::tuple<FourierSeriesTerm, FourierSeriesTerm> multiplication_result;

			PIRANHA_TERM_CTORS(FourierSeriesTerm);

			/// Check if the term is canonical.
			template <class ArgsTuple>
			bool is_canonical(const ArgsTuple &) const 
            {
				return (ancestor::key.sign() > 0);
			}


			// TODO: check if it makes sense to skip the check here and assume canonicalise will be used iff
			// is_canonical has already been tested.
			/// Canonicalise the term.
			template <class ArgsTuple>
			void canonicalise(const ArgsTuple &argsTuple) 
            {
				if (!is_canonical(argsTuple)) 
                {
					invert_trig_args(argsTuple);
				}
			}


			/// Term multiplication.
			/**
			 * NOTE: the result of multiplication here _must_ be canonical.
			 */
			template <class Term1, class Term2, class ArgsTuple>
			static void multiply(const Term1 &t1,
								 const Term2 &t2,
								 multiplication_result &res, const ArgsTuple &argsTuple) 
            {
				// Perform the trigonometric multiplication.
				t1.key.multiply(t2.key, res.template get<0>().key, res.template get<1>().key);
				// Handle coefficient multiplication. Do the first coefficient, then assign the second one.
				// TODO: maybe provide the semantics to coefficients for something like this:
				// cf1.multiply_by_cf(cf2,res.template get<0>().cf,argsTuple),
				// so that we can avoid a copy.
				res.template get<0>().cf = t1.cf;
				res.template get<0>().cf.multBy(t2.cf, argsTuple);
				res.template get<0>().cf.divideBy(2, argsTuple);
				res.template get<1>().cf = res.template get<0>().cf;
				// Now adjust the signs according to werner's formulas.
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

				// Finally, canonicalise the retval terms.
				res.template get<0>().canonicalise(argsTuple);
				res.template get<1>().canonicalise(argsTuple);
			}

		private:

			// Invert the sign of trigonometric multipliers.
			template <class ArgsTuple>
			void invert_trig_args(const ArgsTuple &argsTuple) 
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
