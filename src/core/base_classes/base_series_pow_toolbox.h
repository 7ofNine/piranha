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

#ifndef PIRANHA_BASE_SERIES_POW_TOOLBOX
#define PIRANHA_BASE_SERIES_POW_TOOLBOX

#include <cmath>

#include "../exceptions.h"
#include "../p_assert.h"
#include "../settings_manager.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  template <class Derived>
    struct base_series_pow_toolbox
  {
      template <class ArgsTuple>
        Derived a_pow(const double &x, const ArgsTuple &args_tuple) const throw(unsuitable)
      {
        const int n = (int)nearbyint(x);
        if (std::abs(x - n) <= settings_manager::numerical_zero())
        {
          if (n < 0)
          {
            throw (unsuitable("Cannot raise to negative integer power."));
          }
          else
          {
            Derived retval(natural_pow((size_t)n,args_tuple));
            return retval;
          }
        }
        else
        {
          throw (unsuitable("Cannot raise to real power."));
        }
      }
    private:
      template <class ArgsTuple>
        Derived natural_pow(const size_t &n, const ArgsTuple &args_tuple) const
      {
        typedef typename Derived::term_type term_type;
        typedef typename term_type::cf_type cf_type;
        typedef typename term_type::key_type key_type;
        Derived retval;
        switch (n)
        {
          case 0:
          {
            retval.insert(term_type(cf_type(1,args_tuple),key_type()),args_tuple,derived_const_cast->template nth_index<0>().end());
            break;
          }
          case 1:
          {
            retval.m_container = derived_const_cast->m_container;
            break;
          }
          case 2:
          {
            retval.m_container = derived_const_cast->m_container;
            retval.mult_by(*derived_const_cast,args_tuple);
            break;
          }
          case 3:
          {
            retval.m_container = derived_const_cast->m_container;
            retval.mult_by(*derived_const_cast,args_tuple);
            retval.mult_by(*derived_const_cast,args_tuple);
            break;
          }
          case 4:
          {
            retval.m_container = derived_const_cast->m_container;
            retval.mult_by(*derived_const_cast,args_tuple);
            retval.mult_by(*derived_const_cast,args_tuple);
            retval.mult_by(*derived_const_cast,args_tuple);
            break;
          }
          // Exponentiation by squaring.
          default:
          {
            retval.insert(term_type(cf_type(1,args_tuple),key_type()),args_tuple,derived_const_cast->template nth_index<0>().end());
            Derived tmp(*derived_const_cast);
            size_t i = n;
            while (i)
            {
              if (i & 1)
              {
                retval.mult_by(tmp,args_tuple);
                --i;
              }
              tmp.mult_by(tmp,args_tuple);
              i/=2;
            }
          }
        }
        return retval;
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
