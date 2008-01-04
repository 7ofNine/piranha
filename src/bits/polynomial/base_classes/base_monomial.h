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

#ifndef PIRANHA_BASE_MONOMIAL_H
#define PIRANHA_BASE_MONOMIAL_H

#include <boost/array.hpp> // For ctor from psymbol.

#include "../../psymbol.h"

namespace piranha
{
/// Template parameters for piranha::base_monomial.
#define __BASE_MONOMIAL_TP Cf,Expos
/// Template parameters for piranha::base_monomial (declaration form).
#define __BASE_MONOMIAL_TP_DECL class Cf, class Expos

/// Generica base monomial class.
  template <__BASE_MONOMIAL_TP_DECL>
    class base_monomial
  {
    public:
/// Alias for coefficient type.
      typedef Cf cf_type;
/// Alias for exponents type.
      typedef Expos expos_type;
/// Default ctor.
      explicit base_monomial():m_cf(),m_expos() {}
/// Copy ctor from same or different coeffcient and same exponents.
      template <class Cf2>
        explicit base_monomial(const base_monomial<Cf2,expos_type> &m):m_cf(m.m_cf),m_expos(m.m_expos) {}
/// Ctor from same or different coefficient and same exponents.
      template <class Cf2>
        explicit base_monomial(const Cf2 &cf, const expos_type &e):m_cf(cf),m_expos(e) {}
      explicit base_monomial(const psymbol &);
/// Dtor.
      ~base_monomial() {}
// Getters and setters.
/// Return const reference to coefficient.
      const cf_type *g_cf() const {return &m_cf;}
/// Return const reference to exponents.
      const expos_type *g_expos() const {return &m_expos;}
/// Return mutable reference to coefficient.
      cf_type *g_cf() {return &m_cf;}
/// Return mutable reference to exponents.
      expos_type *g_expos() {return &m_expos;}
// Probing.
    private:
      Cf    m_cf;
      Expo  m_expos;
  };

/// Constructor from psymbol.
  template <__BASE_MONOMIAL_TP_DECL>
    base_monomial<__BASE_MONOMIAL_TP>::base_monomial(const psymbol &p):m_cf(1),m_expos()
  {
    boost::array<int,1> a = { 1 };
    m_expos.assign_int_vector(a);
  }

#undef __BASE_MONOMIAL_TP
#undef __BASE_MONOMIAL_TP_DECL
}

#endif
