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

#include "../base_classes/base_term.h"

namespace piranha
{
  /// Monomial class.
  template <class Cf, class Expo, char Separator, class Allocator>
    class monomial: public base_term<Cf,Expo,Separator,monomial<Cf,Expo,Separator,Allocator>,Allocator>
  {
      /// Alias for the ancestor.
      typedef base_term<Cf,Expo,Separator,monomial<Cf,Expo,Separator,Allocator>,Allocator> ancestor;
      /// Alias for evaluation type.
      typedef typename ancestor::eval_type eval_type;
    public:
      /// Alias for coefficient type.
      typedef Cf cf_type;
      /// Alias for expo type.
      typedef Expo expo_type;
      /// Result of the multiplication of two monomials.
      typedef typename boost::tuple<monomial> multiplication_result;
      /// Default constructor.
      explicit monomial():ancestor() {}
      /// Ctor from string.
      template <class ArgsTuple>
        explicit monomial(const std::string &str, const ArgsTuple &args_tuple):
        ancestor(str,args_tuple)
      {}
      /// Constructor from generic coefficient and fixed exponent part.
      /**
       * Constructs from generic coefficient type.
       */
      template <class Cf2>
        explicit monomial(const Cf2 &c, const expo_type &e):ancestor(cf_type(c),e)
      {}
      /// Generic copy constructor.
      /**
       * Constructs from piranha::monomial with optionally different coefficient type.
       */
      template <class Cf2>
        explicit monomial(const monomial<Cf2,Expo,Separator,Allocator> &term):ancestor(term)
      {}
      /// Monomial multiplication.
      /**
       * NOTE: the result of multiplication here _must_ be canonical.
       */
      template <class Cf2, class ArgsTuple>
        static void multiply(const cf_type &cf1, const expo_type &expo1, const Cf2 &cf2, const expo_type &expo2,
        multiplication_result &res, const ArgsTuple &args_tuple)
      {
        // Perform the multiplication of exponents.
        expo1.multiply(expo2,res.template get<0>().m_key);
        // Handle coefficient multiplication.
        // TODO: maybe provide the semantics to coefficients for something like this:
        // cf1.multiply_by_cf(cf2,res.template get<0>().m_cf,args_tuple),
        // so that we can avoid a copy.
        res.template get<0>().m_cf = cf1;
        res.template get<0>().m_cf.mult_by(cf2,args_tuple);
      }
  };
}

#endif
