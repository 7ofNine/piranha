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

namespace piranha
{
  /// Natural power for series.
  /**
   * Calculated through exponentiation by squaring.
   */
  template <class T, class ArgsTuple>
    T natural_power(const T &x, const size_t &n, const ArgsTuple &args_tuple)
  {
    T retval;
    switch (n)
    {
      case 0:
      {
        retval = T(1,args_tuple);
        break;
      }
      case 1:
      {
        retval = x;
        break;
      }
      case 2:
      {
        retval = x;
        retval.mult_by(x,args_tuple);
        break;
      }
      case 3:
      {
        retval = x;
        retval.mult_by(x,args_tuple);
        retval.mult_by(x,args_tuple);
        break;
      }
      case 4:
      {
        retval = x;
        retval.mult_by(x,args_tuple);
        retval.mult_by(x,args_tuple);
        retval.mult_by(x,args_tuple);
        break;
      }
      default:
      {
        retval = T(1,args_tuple);
        // Use scoping here to have tmp destroyed when it is not needed anymore.
        {
          T tmp(x);
          size_t i = n;
          while (i)
          {
            if (i & 1)
            {
              retval.mult_by(tmp,args_tuple);
              --i;
            }
            i/=2;
            if (i != 0)
            {
              tmp.mult_by(tmp,args_tuple);
            }
          }
        }
      }
    }
    return retval;
  }
}

#endif
