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
#include <boost/type_traits/is_same.hpp>
#include <cstddef>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../settings.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Binomial exponentiation toolbox.
	/**
	 * Overrides BaseSeries::real_power, BaseSeries::negative_integer_power, BaseSeries::rational_power
	 * and reimplements them using binomial expansion.
	 */
	template <class Derived>
	class binomial_exponentiation
	{
		public:

			/// Real power.
			template <class ArgsTuple>
			Derived real_power(const double &y, const ArgsTuple &argsTuple) const
			{
				return generic_binomial_power(
					derived_const_cast->template get_sorted_series<Derived>(argsTuple),y,argsTuple);
			}


			/// Negative integer power.
			template <class ArgsTuple>
			Derived negative_integer_power(const int &y, const ArgsTuple &argsTuple) const
			{
				return generic_binomial_power(
					derived_const_cast->template get_sorted_series<Derived>(argsTuple),y, argsTuple);
			}


			/// Rational power.
			template <class ArgsTuple>
			Derived rational_power(const mp_rational &q, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(q != 0 && q != 1);
				return generic_binomial_power(
					derived_const_cast->template get_sorted_series<Derived>(argsTuple), q, argsTuple);
			}

		private:

			template <class Term, class Number, class ArgsTuple>
			static Derived generic_binomial_power(const std::vector<Term const *> &v, const Number &y, const ArgsTuple &argsTuple)
			{
				typedef typename Derived::TermType term_type;
				// Here we know that the cases of empty series and natural power have already
				// been taken care of in BaseSeries::base_pow.
				PIRANHA_ASSERT(v.size() >= 1);

				term_type A(*v[0]);
				// This is X, i.e., the original series without the leading term, which will then be divided by A.
				Derived XoverA;
				const std::size_t size = v.size();
				for (std::size_t i = 1; i < size; ++i) 
                {
					XoverA.insert(term_type(*v[i]),argsTuple);
				}

				// Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
				term_type tmp_term(A.cf.pow(-1, argsTuple), A.key.pow(-1, argsTuple));
				Derived Ainv;
				Ainv.insert(tmp_term, argsTuple);
				// Now let's compute X/A.
				XoverA.base_mult_by(Ainv, argsTuple);
				// Get the expansion limit from the truncator.
				std::size_t n;
				try {
					n = XoverA.psi_(0, 1, argsTuple);
				} catch (const value_error &ve) 
                {
					PIRANHA_THROW(value_error,std::string("series is unsuitable for exponentiation through binomial expansion."
						"\nThe reported error is: ")
						+ ve.what());
				}

				return binomial_expansion(A, XoverA, y, n, argsTuple);
			}


			template <class Term, class Number, class ArgsTuple>
			static Derived binomial_expansion(const Term &A, const Derived &XoverA,
				const Number &y, const std::size_t &n, const ArgsTuple &argsTuple)
			{
				typedef typename Derived::TermType term_type;

				PIRANHA_STATIC_CHECK((boost::is_same<Term, typename Derived::TermType>::value),
					                 "Term type mismatch in binomial expansion.");

				// Start the binomial expansion.
				term_type tmp_term;
				// Calculate A**y. See if we can raise to real power the coefficient and the key.
				// Exceptions will be thrown in case of problems.
				tmp_term.cf  = A.cf.pow(y, argsTuple);
				tmp_term.key = A.key.pow(y, argsTuple);
				Derived Apowy;
				Apowy.insert(tmp_term, argsTuple);
				// Let's proceed now to the bulk of the binomial expansion. Luckily we can compute the needed generalised
				// binomial coefficient incrementally at every step. We start with 1.
				Derived retval;
				Derived tmp;
				tmp.baseAdd(1, argsTuple);
				retval.baseAdd(tmp, argsTuple);
				for (std::size_t i = 1; i < n; ++i) 
                {
					tmp.base_mult_by(y - (double)i + 1, argsTuple);
					tmp.base_divide_by(boost::numeric_cast<int>(i), argsTuple);
					tmp.base_mult_by(XoverA, argsTuple);
					retval.baseAdd(tmp, argsTuple);
				}
				// Finally, multiply the result of the summation by A**y.
				retval.base_mult_by(Apowy, argsTuple);
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
