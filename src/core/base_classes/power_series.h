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

#ifndef PIRANHA_POWER_SERIES_H
#define PIRANHA_POWER_SERIES_H

#include <boost/static_assert.hpp>

#include "../integer_typedefs.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  /// Power series toolbox.
  /**
  * This toolbox assumes that the terms of the series are sorted in ascending total degree.
  */
  template <int ExpoPosition, class Derived>
    class power_series
  {
      BOOST_STATIC_ASSERT(ExpoPosition >= 0);
    public:
      /// Get the degree of the power series.
      max_fast_int degree() const
      {
        if (derived_const_cast->template nth_index<0>().empty())
        {
          return 0;
        }
        typename Derived::const_sorted_iterator it = derived_const_cast->template nth_index<0>().end();
        --it;
        return (it->m_cf.degree()+it->m_key.degree());
      }
      /// Get the minimum degree of the power series.
      max_fast_int min_degree() const
      {
        if (derived_const_cast->template nth_index<0>().empty())
        {
          return 0;
        }
        return derived_const_cast->template nth_index<0>().begin()->m_cf.min_degree()+
          derived_const_cast->template nth_index<0>().begin()->m_key.min_degree();
      }
      static const int expo_position = ExpoPosition;
      void upload_min_exponents(std::vector<max_fast_int> &v) const
      {
        p_assert(!derived_const_cast->empty());
        typedef typename Derived::const_sorted_iterator const_sorted_iterator;
        const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
        const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin();
        it->m_cf.upload_min_exponents(v);
        it->m_key.upload_min_exponents(v);
        for (; it != it_f; ++it)
        {
          it->m_cf.test_min_exponents(v);
          it->m_key.test_min_exponents(v);
        }
      }
      /// Get a vector containing the minimum exponents of the series.
      template <class ArgsTuple>
        std::vector<max_fast_int> min_exponents(const ArgsTuple &args_tuple) const
      {
        std::vector<max_fast_int> retval(args_tuple.template get<ExpoPosition>().size());
        upload_min_exponents(retval);
        return retval;
      }
      void test_min_exponents(std::vector<max_fast_int> &v) const
      {
        typedef typename Derived::const_sorted_iterator const_sorted_iterator;
        const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
        for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it)
        {
          it->m_cf.test_min_exponents(v);
          it->m_key.test_min_exponents(v);
        }
      }
      // Return true if the minimum exponents are smaller than those specified in the limits vector.
      template <class ArgsTuple>
        bool test_expo_limits(const std::vector<std::pair<size_t,max_fast_int> > &v, const ArgsTuple &args_tuple) const
      {
        const std::vector<max_fast_int> min_expo(min_exponents(args_tuple));
        const size_t size = v.size();
        for (size_t i = 0; i < size; ++i)
        {
          p_assert(v[i].first < min_expo.size());
          if (v[i].second < min_expo[v[i].first])
          {
            return false;
          }
        }
        return true;
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
