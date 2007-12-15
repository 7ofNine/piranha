/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_NORM_TRUNCATION_TOOLBOX_H
#define PIRANHA_NORM_TRUNCATION_TOOLBOX_H

#include <iostream>

namespace piranha
{
  class norm_truncation_toolbox
  {
    public:
      static void set_truncation(const double &x)
      {
        if (x < 0)
        {
          std::cout << "Please provide a non-negative value." << std::endl;
        }
        else
        {
          trunc_level = x;
        }
      }
      static const double &get_truncation()
      {
        return trunc_level;
      }
    private:
      static double trunc_level;
  };
}

#endif
