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

#include <boost/type_traits/is_same.hpp>
#include <string>

#include "../exceptions.h"
#include "../p_assert.h"
#include "../settings.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class binomial_exponentiation_toolbox
	{
		protected:
			/// Real power.
			/**
			* This method is written to work in conjunction with base_series::b_pow.
			*/
			template <class ArgsTuple>
			Derived real_power(const double &y, const ArgsTuple &args_tuple) const {
				return generic_power(y, args_tuple);
			}
			template <class ArgsTuple>
			Derived negative_integer_power(const max_fast_int &y, const ArgsTuple &args_tuple) const {
				return generic_power(y, args_tuple);
			}
		private:
			template <class Number, class ArgsTuple>
			Derived generic_power(const Number &y, const ArgsTuple &args_tuple) const {
				typedef typename Derived::term_type term_type;
				// Here we know that the cases of single term, empty series and natural power have already been taken care of
				// in base_series::b_pow.
				p_assert(derived_const_cast->template nth_index<0>().size() > 1);
				term_type A(*derived_const_cast->template nth_index<0>().begin());
				// This is X, i.e., the original series without the leading term, which will then be divided by A.
				Derived XoverA(*derived_const_cast);
				XoverA.template term_erase<0>(args_tuple, XoverA.template nth_index<0>().begin());
				// Now let's try to calculate 1/A. There will be exceptions thrown if we cannot do that.
				term_type tmp_term(A.m_cf.pow((max_fast_int)(-1), args_tuple), A.m_key.pow(max_fast_int(-1), args_tuple));
				Derived Ainv;
				Ainv.insert(tmp_term, args_tuple, Ainv.template nth_index<0>().end());
				// Now let's compute X/A.
				XoverA.mult_by(Ainv, args_tuple);
				// Get the expansion limit from the truncator.
				size_t n;
				try {
					n = Derived::multiplier_type::truncator_type::power_series_limit(XoverA, args_tuple);
				} catch (const unsuitable &u) {
					throw unsuitable(std::string("Series is unsuitable for real exponentiation through binomial expansion."
						"\nThe reported error is: ")
									 + u.what());
				}
				return binomial_expansion(A, XoverA, y, n, args_tuple);
			}
			template <class Term, class Number, class ArgsTuple>
			static Derived binomial_expansion(const Term &A, const Derived &XoverA,
									const Number &y, const size_t &n, const ArgsTuple &args_tuple)
			{
				typedef typename Derived::term_type term_type;
				BOOST_STATIC_ASSERT((boost::is_same<term_type, typename Derived::term_type>::value));
				// Start the binomial expansion.
				term_type tmp_term;
				// Calculate A**y. See if we can raise to real power the coefficient and the key.
				// Exceptions will be thrown in case of problems.
				tmp_term.m_cf = A.m_cf.pow(y, args_tuple);
				tmp_term.m_key = A.m_key.pow(y, args_tuple);
				Derived Apowy;
				Apowy.insert(tmp_term, args_tuple, Apowy.template nth_index<0>().end());
				// Let's proceed now to the bulk of the binomial expansion. Luckily we can compute the needed generalised
				// binomial coefficient incrementally at every step. We start with 1.
				Derived retval;
				Derived tmp((max_fast_int)1, args_tuple);
				retval.add(tmp, args_tuple);
				for (size_t i = 1; i <= n; ++i) {
					tmp.mult_by(y-(max_fast_int)(i)+(max_fast_int)1, args_tuple);
					tmp.divide_by((max_fast_int)i, args_tuple);
					tmp.mult_by(XoverA, args_tuple);
					retval.add(tmp, args_tuple);
				}
				// Finally, multiply the result of the summation by A**y.
				retval.mult_by(Apowy, args_tuple);
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
