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

#ifndef PIRANHA_CODED_SERIES_MULTIPLIER_H
#define PIRANHA_CODED_SERIES_MULTIPLIER_H

#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <valarray>

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  template <class Derived>
    class coded_series_multiplier
  {
    protected:
      coded_series_multiplier():
        m_cr_is_viable(false),
        m_size(derived_const_cast->m_args_tuple.template get<Derived::key_type::position>().size()),
        m_min_max1(m_size),m_min_max2(m_size),m_res_min_max(m_size),m_fast_res_min_max(m_size)
      {}
    protected:
      // Is coded representation viable?
      bool                                                  m_cr_is_viable;
      // Size of the coding vector, min_max vectors, etc.
      const size_t                                          m_size;
      // Vectors of minimum and maximum value pairs for the series being multiplied.
      std::valarray<std::pair<max_fast_int,max_fast_int> >  m_min_max1;
      std::valarray<std::pair<max_fast_int,max_fast_int> >  m_min_max2;
      // Vector of minimum and maximum value pairs for the resulting series.
      // GMP is used to avoid trespassing the range limits of max_fast_int.
      std::valarray<std::pair<mpz_class,mpz_class> >        m_res_min_max;
      // Version of the above downcast to fast integer type.
      std::valarray<std::pair<max_fast_int,max_fast_int> >  m_fast_res_min_max;
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
