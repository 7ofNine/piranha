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

#ifndef PIRANHA_NORM_TRUNCATABLE_SERIES_H
#define PIRANHA_NORM_TRUNCATABLE_SERIES_H

#include <cmath>
#include <iostream>

#include "../config.h"

namespace piranha
{
  class norm_truncatable_series
  {
    public:
      static const double &get_truncation()
      {
        return m_truncation_level;
      }
      static void set_truncation(const int &n)
      {
        if (n < 0)
        {
          std::cout << "Please insert a non-negative integer." << std::endl;
        }
        else if (n == 0)
        {
          m_truncation_level = 0;
        }
        else
        {
          m_truncation_level = std::pow(10,-n);
        }
      }
    private:
      static double m_truncation_level;
  };
}

#endif
