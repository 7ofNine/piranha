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

#ifndef PIRANHA_BASE_SERIES_SPECIAL_FUNCTIONS_H
#define PIRANHA_BASE_SERIES_SPECIAL_FUNCTIONS_H


#include <cstddef>
#include <string>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#include "../exceptions.h"
#include "../math.h"
#include "../mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class BaseSeriesSpecialFunctions
	{
		//protected:
		public:
            // hypergeometric
			template <class ArgsTuple>
			Derived baseHyperF(std::vector<mp_rational> const &aList, std::vector<mp_rational> const &bList,  int const iterationLimit, ArgsTuple const &argsTuple) const
			{
				Derived retval;
				// If one of the elements in bList is a non-positive integer, a division by zero will occur.
				hyperFbCheck(bList);
				// Cache values.
				const std::size_t aSize = aList.size();
                const std::size_t bSize = bList.size();
				// HyperF of a null series will always be equal to 1.
				if (derived_const_cast->empty()) 
				{
					retval.baseAdd(1, argsTuple);
					return retval;
				}

				// If one of the elements in a_list is zero, hyperF will always be 1.
				for (std::size_t i = 0; i < aSize; ++i) 
				{
					if (aList[i] == 0) 
					{
						retval.baseAdd(1, argsTuple);
						return retval;
					}
				}

				// Establish series limit, using the following strategy...
				int iter;
				try {
					// First we try to respect the given iterationLimit, if any,
					// otherwise we call the truncator to see what it has to say.
					iter = (iterationLimit >= 0) ? iterationLimit : derived_const_cast->psIterations(0, 1, argsTuple);

				} catch (const value_error &) 
				{
					// If the truncator fails to establish a limit, we go to look into aList and bList
					// whether there are negative integer numbers.
					mp_rational const *minInt = 0;
					for (std::size_t i = 0; i < aSize; ++i) 
					{
						if (aList[i] < 0 && aList[i].get_den() == 1) 
						{
							if (!minInt || aList[i] < *minInt) 
							{
								minInt = &aList[i];
							}
						}
					}

					if (!minInt) 
					{
						PIRANHA_THROW(value_error, "Could not establish a value for the series limit in hyperF: "
							"no explicit limit was provided, the truncator failed to establish a limit and "
							"no negative integer numbers were found in aList");
					}

					PIRANHA_ASSERT(*minInt < 0);

					iter = -(((*minInt) - 1).to_int());
				}

				const std::size_t maxIter = iter;
				// If no iterations are needed, return an empty series.
				if (!maxIter) 
				{
					return retval;
				}

				retval.baseAdd(1, argsTuple);
				Derived tmp;
				tmp.baseAdd(1, argsTuple);
				for (std::size_t i = 1; i < maxIter; ++i) 
				{
					for (std::size_t j = 0; j < aSize; ++j) 
					{
						tmp.baseMultBy(aList[j] + (int)(i - 1), argsTuple);
					}

					for (std::size_t j = 0; j < bSize; ++j) 
					{
						tmp.baseDivideBy(bList[j] + (int)(i - 1), argsTuple);
					}

					tmp.baseDivideBy(boost::numeric_cast<int>(i), argsTuple);
					tmp.baseMultBy(*derived_const_cast, argsTuple);
					retval.baseAdd(tmp, argsTuple);
				}

				return retval;
			}

            // n-th derivateive of hypergeoemtric series
			template <class ArgsTuple>
			Derived baseDHyperF(const int n, const std::vector<mp_rational> &aList, const std::vector<mp_rational> &bList, const int &iterationLimit, const ArgsTuple &argsTuple) const
			{
				if (n < 0) 
				{
					PIRANHA_THROW(value_error, "please enter a non-negative order of differentiation");
				}

				// If one of the elements in b_list is a non-positive integer, a division by zero will occur.
				hyperFbCheck(bList);
				// Compute the outside factor.
				mp_rational factor(1);
				const std::size_t aSize = aList.size();
                const std::size_t bSize = bList.size();
				for (std::size_t i = 0; i < aSize; ++i) 
				{
					factor *= aList[i].r_factorial(n);
				}

				for (std::size_t i = 0; i < bSize; ++i) 
				{
					factor /= bList[i].r_factorial(n);
				}

				// Create the new a_list and b_list.
				std::vector<mp_rational> aListNew(aList);
                std::vector<mp_rational> bListNew(bList);
				for (std::size_t i = 0; i < aSize; ++i) 
				{
					aListNew[i] += n;
				}

				for (std::size_t i = 0; i < bSize; ++i)
                {
					bListNew[i] += n;
				}

				Derived retval(baseHyperF(aListNew, bListNew, iterationLimit, argsTuple));
				retval.baseMultBy(factor, argsTuple);
				return retval;
			}


			/// Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived baseBesselJ(const int order, const ArgsTuple &argsTuple) const 
			{
				Derived retval;
				// Dispatch besselJ to coefficient if series consists of a single coefficient only.
				if (derived_const_cast->isSingleCf()) 
				{
					typedef typename Derived::TermType termType;
					retval.insert(termType(derived_const_cast->begin()->cf.besselJ(order, argsTuple), typename termType::KeyType()), argsTuple);
					return retval;
				}

				// Take care of negative order.
				const int orderNormalized = (order >= 0) ? order : - order;
				// Calculate this/2.
				Derived x_2(*derived_const_cast);
				x_2.baseDivideBy(2, argsTuple);

				retval = Derived(x_2).baseMultBy(-1, argsTuple).baseMultBy(x_2, argsTuple)
					                 .baseHyperF(std::vector<mp_rational>(), std::vector<mp_rational>((std::size_t)1, mp_rational(orderNormalized) + 1), -1, argsTuple);
				retval.baseMultBy(x_2.basePow(orderNormalized, argsTuple), argsTuple);
				retval.baseDivideBy(mp_integer(orderNormalized).factorial(), argsTuple);

				if (order < 0) 
				{
					retval.baseMultBy(cs_phase(order), argsTuple); //(-1)^n
				}

				return retval;
			}


			/// Partial derivative with respect to the argument of Bessel function of the first kind of integer order.
			template <class ArgsTuple>
			Derived baseDBesselJ(const int order, const ArgsTuple &argsTuple) const
			{
				// NOTE: we don't user hyperF here because direct series expansion is faster.
				const int orderNormalized = (order >= 0) ? order : - order;
				Derived retval;

				if (!orderNormalized) 
				{
					// Exploit the fact that for order == 0, derivative is -J_1.
					retval = baseBesselJ(1, argsTuple);
					retval.baseMultBy(-1, argsTuple);
					return retval;
				}

				// Get the expansion limit from the truncator.
				int limit = 0;
				try {
					limit = derived_const_cast->psIterations(orderNormalized - 1, 2, argsTuple);
				} catch (const value_error &ve) 
				{
					PIRANHA_THROW(value_error,std::string("series is unsuitable as argument of the derivative of "
						"Bessel function of the first kind.\nThe reported error is: ")
						+ ve.what());
				}

				if (!limit) 
				{
					return retval;
				}

				// Now we buid the starting point of the power series expansion.
				retval = *derived_const_cast;
				retval.baseDivideBy(2, argsTuple);
				
                // This will be used later.
				Derived square_x2(retval);
				square_x2.baseMultBy(square_x2, argsTuple);
				retval = retval.basePow(orderNormalized - 1, argsTuple);
				retval.baseDivideBy(mp_integer(orderNormalized).factorial(), argsTuple);
				retval.baseMultBy(orderNormalized, argsTuple);
				retval.baseDivideBy(2, argsTuple);
				
                // Now let's proceed to the bulk of the power series expansion.
				Derived tmp(retval);

				for (int i = 1; i < (int)limit; ++i) 
				{
					tmp.baseMultBy((-1) * (mp_integer(orderNormalized) + 2 * mp_integer(i)), argsTuple);
					tmp.baseDivideBy((mp_integer(i) * (mp_integer(i) + orderNormalized)) * (mp_integer(orderNormalized) + 2 * (mp_integer(i) - 1)), argsTuple);
					tmp.baseMultBy(square_x2, argsTuple);
					retval.baseAdd(tmp, argsTuple);
				}

				return retval;
			}


			/// Bessel function of the first kind of integer order divided by its argument**m.
			template <class ArgsTuple>
			Derived baseBesselJDivm(const int order, const int m, const ArgsTuple &argsTuple) const
			{
				// Take care of negative order.
				const int orderNormalized = (order >= 0) ? order : - order;
				if (orderNormalized < m) 
				{
					PIRANHA_THROW(zero_division_error, "absolute value of order must not be smaller than m");
				}

				Derived retval;
				// Calculate this/2.
				Derived x_2(*derived_const_cast);
				x_2.baseDivideBy(2, argsTuple);

				retval = Derived(x_2).baseMultBy(-1, argsTuple).baseMultBy(x_2, argsTuple)
					                 .baseHyperF(std::vector<mp_rational>(), std::vector<mp_rational>((std::size_t)1, mp_rational(orderNormalized) + 1), -1, argsTuple);
				retval.baseMultBy(derived_const_cast->basePow(orderNormalized - m, argsTuple), argsTuple);
				retval.baseDivideBy(mp_integer(orderNormalized).factorial() * mp_integer(2).pow(orderNormalized), argsTuple);

				if (order < 0) 
				{
					retval.baseMultBy(cs_phase(order), argsTuple);
				}

				return retval;
			}

		private:

            //check the hypergeometric b-list for negative numbers or 0 . Not good to have it
			static void hyperFbCheck(const std::vector<mp_rational> &bList)
			{
				const std::size_t bSize = bList.size();
				for (std::size_t i = 0; i < bSize; ++i) 
				{
					if (bList[i] <= 0 && bList[i].get_den() == 1) 
					{
						PIRANHA_THROW(zero_division_error, "bList in hyperF contains a non-positive integer");
					}
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
