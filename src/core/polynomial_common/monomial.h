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

#ifndef PIRANHA_MONOMIAL_H
#define PIRANHA_MONOMIAL_H

#include <boost/tuple/tuple.hpp>

#include "../base_classes/base_term.h"

namespace piranha
{
	/// Monomial class.
	template <class Cf, class Key, char Separator, class Allocator>
	class Monomial: public BaseTerm<Cf, Key, Separator, Allocator, Monomial<Cf, Key, Separator, Allocator> >
	{
			// Alias for the ancestor.
			typedef BaseTerm<Cf, Key, Separator, Allocator, Monomial> ancestor;

		public:

			/// Alias for coefficient type.
			typedef Cf CfType;
			/// Alias for expo type.
			typedef Key KeyType;
			/// Result of the multiplication of two monomials.
			typedef typename boost::tuple<Monomial> multiplication_result;

			PIRANHA_TERM_CTORS(Monomial);

			/// Monomial multiplication.
			/**
			 * NOTE: the result of multiplication here _must_ be canonical.
			 */
			template <class Term1, class Term2, class ArgsTuple>
			static void multiply(const Term1 &m1, const Term2 &m2, multiplication_result &res, const ArgsTuple &argsTuple) 
            {
				// Perform the multiplication of exponents.
				m1.key.multiply(m2.key, res.template get<0>().key);
				// Handle coefficient multiplication.
				// TODO: maybe provide the semantics to coefficients for something like this:
				// cf1.multiply_by_cf(cf2,res.template get<0>().cf,argsTuple),
				// so that we can avoid a copy.
				res.template get<0>().cf = m1.cf;
				res.template get<0>().cf.multBy(m2.cf, argsTuple);
			}
	};
}

#endif
