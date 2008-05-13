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

#ifndef PIRANHA_SERIES_MATH_H
#define PIRANHA_SERIES_MATH_H

/*! \file series_math.h
    \brief Mathematics for series.

    The functions defined here require series and arguments tuples as parameters.
*/

#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.

#include "../exceptions.h"
#include "../integer_typedefs.h"

namespace piranha
{
	/// Natural power for series.
	/**
	 * Calculated through exponentiation by squaring.
	 */
	template <class Series, class ArgsTuple>
	Series natural_power(const Series &x, const size_t &n, const ArgsTuple &args_tuple)
	{
		Series retval;
		switch (n) {
		case 0: {
			retval = Series((max_fast_int)1, args_tuple);
			break;
		}
		case 1: {
			retval = x;
			break;
		}
		case 2: {
			retval = x;
			retval.mult_by(x, args_tuple);
			break;
		}
		case 3: {
			retval = x;
			retval.mult_by(x, args_tuple);
			retval.mult_by(x, args_tuple);
			break;
		}
		case 4: {
			retval = x;
			retval.mult_by(x, args_tuple);
			retval.mult_by(x, args_tuple);
			retval.mult_by(x, args_tuple);
			break;
		}
		default: {
			retval = Series((max_fast_int)1, args_tuple);
			// Use scoping here to have tmp destroyed when it is not needed anymore.
			{
				Series tmp(x);
				size_t i = n;
				while (i) {
					if (i & 1) {
						retval.mult_by(tmp, args_tuple);
						--i;
					}
					i /= 2;
					if (i != 0) {
						tmp.mult_by(tmp, args_tuple);
					}
				}
			}
		}
		}
		return retval;
	}

	template <class Series, class ArgsTuple>
	Series binomial_expansion(const typename Series::term_type &A, const Series &XoverA,
							  const double &y, const size_t &n, const ArgsTuple &args_tuple)
	{
		typedef typename Series::term_type term_type;
		// Start the binomial expansion.
		term_type tmp_term;
		// Calculate A**y. See if we can raise to real power the coefficient and the key.
		// Exceptions will be thrown in case of problems.
		tmp_term.m_cf = A.m_cf.pow(y, args_tuple);
		tmp_term.m_key = A.m_key.pow(y, args_tuple);
		Series Apowy;
		Apowy.insert(tmp_term, args_tuple, Apowy.template nth_index<0>().end());
		// Let's proceed now to the bulk of the binomial expansion. Luckily we can compute the needed generalised
		// binomial coefficient incrementally at every step. We start with 1.
		Series retval;
		Series tmp((max_fast_int)1, args_tuple);
		retval.add(tmp, args_tuple);
		mpq_class mpq_y;
		for (size_t i = 1; i <= n; ++i) {
			// This hack is to make sure we transform y into a fraction as early as possible,
			// and to do the calculations in rational form. This way we should
			// minimize the chance of losing precision due to conversion to/from floating point.
			mpq_y = y;
			mpq_y -= i;
			mpq_y += 1;
			tmp.mult_by(mpq_y.get_d(), args_tuple);
			tmp.divide_by((max_fast_int)i, args_tuple);
			tmp.mult_by(XoverA, args_tuple);
			retval.add(tmp, args_tuple);
		}
		// Finally, multiply the result of the summation by A**y.
		retval.mult_by(Apowy, args_tuple);
		return retval;
	}
}

#endif
