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

#ifndef PIRANHA_SERIES_POW_TOOLBOX
#define PIRANHA_SERIES_POW_TOOLBOX

#include <cmath>

#include "../exceptions.h"
#include "../settings_manager.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  template <class Derived>
    struct series_pow_toolbox
  {
    public:
      template <class ArgsTuple>
        Derived pow(const double &x, const ArgsTuple &args_tuple) const throw(unsuitable)
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
            Derived retval = natural_pow((size_t)n,args_tuple);
            return retval;
          }
        }
        else
        {
          throw (unsuitable("Cannot raise to real power."));
        }
      }
      Derived pow(const double &x) const throw(unsuitable)
      {
        Derived retval = pow(x,derived_const_cast->m_arguments);
        return retval;
      }
    private:
      template <class ArgsTuple>
        Derived natural_pow(const size_t &n, const ArgsTuple &args_tuple) const
      {
        switch (n)
        {
          case 0:
          {
            Derived retval(1,args_tuple);
            return retval;
          }
          case 1:
          {
            Derived retval(*derived_const_cast);
            return retval;
          }
          default:
          {
            Derived retval(*derived_const_cast);
            for (size_t i = 1; i < n; ++i)
            {
              retval.mult_by((*derived_const_cast),args_tuple);
            }
            return retval;
          }
        }
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
