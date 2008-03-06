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
      /// Default constructor.
      explicit monomial():ancestor::base_term() {}
      /// Ctor from string.
      template <class ArgsTuple>
        explicit monomial(const std::string &str, const ArgsTuple &args_tuple):
        ancestor::base_term(str,args_tuple)
      {}
      /// Constructor from generic coefficient and fixed exponent part.
      /**
       * Constructs from generic coefficient type.
       */
      template <class Cf2>
        explicit monomial(const Cf2 &c, const expo_type &e):ancestor(cf_type(c),t)
      {}
      /// Generic copy constructor.
      /**
       * Constructs from piranha::monomial with optionally different coefficient type.
       */
      template <class Cf2>
        explicit monomial(const monomial<Cf2,Expo,Separator,Allocator> &term):ancestor(term)
      {}
  };
}

#endif
