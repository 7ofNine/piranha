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

#ifndef PIRANHA_BINOMIAL_EXPONENTIATION_TOOLBOX_H
#define PIRANHA_BINOMIAL_EXPONENTIATION_TOOLBOX_H

#include <boost/numeric/conversion/cast.hpp>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../settings.h"

#include <cstddef>
#include <string>
#include <vector>
#include <type_traits>

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Binomial exponentiation toolbox.
	/**
	 * Overrides BaseSeries::realPower, BaseSeries::negativeIntegerPower, BaseSeries::rationalPower
	 * and reimplements them using binomial expansion.
     * 
	 */
    //
    // requires the SeriesMultiplication and truncator for the DerivedSeries
    //

	template <class Derived>
	class BinomialExponentiation
	{
		public:

			/// Real power.
			template <class ArgsTuple>
			Derived realPower(const double y, const ArgsTuple &argsTuple) const
			{
				return genericBinomialPower(derived_const_cast->template get_sorted_series<Derived>(argsTuple), y, argsTuple);
			}


			/// Negative integer power.
			template <class ArgsTuple>
			Derived negativeIntegerPower(const int y, const ArgsTuple &argsTuple) const
			{
				return genericBinomialPower(derived_const_cast->template get_sorted_series<Derived>(argsTuple), y, argsTuple);
			}


			/// Rational power.
			template <class ArgsTuple>
			Derived rationalPower(const mp_rational &q, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(q != 0 && q != 1);

				return genericBinomialPower(derived_const_cast->template get_sorted_series<Derived>(argsTuple), q, argsTuple);
			}

		private:

			template <class Term, class Number, class ArgsTuple>
			static Derived genericBinomialPower(const std::vector<Term const *> &v, const Number &y, const ArgsTuple &argsTuple)
			{
				typedef typename Derived::TermType TermType;
				// Here we know that the cases of empty series and natural power have already
				// been taken care of in BaseSeries::basePow.
				PIRANHA_ASSERT(v.size() >= 1);

				TermType a(*v[0]);
				// This is x, i.e., the original series without the leading term, which will then be divided by A.
				Derived x_a;
				const std::size_t size = v.size();
				for (std::size_t i = 1; i < size; ++i) 
                {
					x_a.insert(TermType(*v[i]), argsTuple);
				}

				// Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
				TermType tmpTerm(a.cf.pow(-1, argsTuple), a.key.pow(-1, argsTuple));
				Derived aInv;
				aInv.insert(tmpTerm, argsTuple);
				// Now let's compute X/A.
				x_a.baseMultBy(aInv, argsTuple);
				// Get the expansion limit from the truncator.
				int n;
				try {
					n = x_a.psIterations(0, 1, argsTuple);

				} catch (const value_error &ve) 
                {
					PIRANHA_THROW(value_error, std::string("Series is unsuitable for exponentiation through binomial expansion."
						"\nThe reported error is: ")
						+ ve.what());
				}

				return binomialExpansion(a, x_a, y, n, argsTuple);
			}


			template <class Term, class Number, class ArgsTuple>
			static Derived binomialExpansion(const Term &a, const Derived &x_a, const Number &y, const std::size_t n, const ArgsTuple &argsTuple)
			{
				typedef typename Derived::TermType TermType;

                static_assert((std::is_same_v<Term, typename Derived::TermType>), "Term type mismatch in binomial expansion.");

				// Start the binomial expansion.
				TermType tmpTerm;
				// Calculate a**y. See if we can raise to real power the coefficient and the key.
				// Exceptions will be thrown in case of problems.
				tmpTerm.cf  = a.cf.pow(y, argsTuple);
				tmpTerm.key = a.key.pow(y, argsTuple);
				Derived ay;
				ay.insert(tmpTerm, argsTuple);
				// Let's proceed now to the bulk of the binomial expansion. Luckily we can compute the needed generalised
				// binomial coefficient incrementally at every step. We start with 1.
				Derived retval;
				Derived tmp;
				tmp.baseAdd(1, argsTuple);
				retval.baseAdd(tmp, argsTuple);

				for (std::size_t i = 1; i < n; ++i) 
                {
					tmp.baseMultBy(y - (double)i + 1, argsTuple);
					tmp.baseDivideBy(boost::numeric_cast<int>(i), argsTuple);
					tmp.baseMultBy(x_a, argsTuple);
					retval.baseAdd(tmp, argsTuple);
				}

				// Finally, multiply the result of the summation by A**y.
				retval.baseMultBy(ay, argsTuple);

				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
