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

#include "../p_assert.h"

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
        m_min_max1(m_size),m_min_max2(m_size),m_res_min_max(m_size),m_fast_res_min_max(m_size),
        // Coding vector is larger to accomodate extra element at the end.
        m_coding_vector(m_size+1)
      {}
      void find_input_min_max()
      {
        typedef typename Derived::iterator1 iterator1;
        typedef typename Derived::iterator2 iterator2;
        const iterator1 it_f1 = derived_cast->m_s1.template nth_index<0>().end();
        const iterator2 it_f2 = derived_cast->m_s2.template nth_index<0>().end();
        iterator1 it1 = derived_cast->m_s1.template nth_index<0>().begin();
        iterator2 it2 = derived_cast->m_s2.template nth_index<0>().begin();
        // Fill first minmax vector. This works because at this point we are sure both series have
        // at least one term. Assert it, just to make sure.
        p_assert(!derived_cast->m_s1.template nth_index<0>().empty() and !derived_cast->m_s2.template nth_index<0>().empty());
        it1->m_key.template update_limits<true>(m_min_max1);
        it2->m_key.template update_limits<true>(m_min_max2);
        // Move to the second terms and cycle on all remaining terms.
        ++it1;
        ++it2;
        for (; it1 != it_f1; ++it1)
        {
          it1->m_key.template update_limits<false>(m_min_max1);
        }
        for (; it2 != it_f2; ++it2)
        {
          it2->m_key.template update_limits<false>(m_min_max2);
        }
std::cout << "Limits are:\n";
for (size_t i = 0; i < m_min_max1.size(); ++i)
{
  std::cout << m_min_max1[i].first << ',' << m_min_max1[i].second << '\n';
}
std::cout << "and:\n";
for (size_t i = 0; i < m_min_max2.size(); ++i)
{
  std::cout << m_min_max2[i].first << ',' << m_min_max2[i].second << '\n';
}
      }
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
      // Coding vector.
      std::valarray<max_fast_int>                           m_coding_vector;
      // Mininum and maximum values of codes.
      max_fast_int                                          m_h_min;
      max_fast_int                                          m_h_max;
      // Coded keys.
      std::valarray<max_fast_int>                           m_ckeys1;
      std::valarray<max_fast_int>                           m_ckeys2;
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
