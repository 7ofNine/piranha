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

#include "symbolic_truncator.h"

namespace piranha
{
  vec_expo_limit symbolic_truncator::vel_;

/// Check whether a limit has already been set for a symbol.
/**
 * Positive result means found, negative means not found.
 */
  int symbolic_truncator::find_expo_limit(psym_p it)
  {
    const size_t w=vel_.size();
    for (size_t j=0;j<w;++j)
      {
        if (vel_[j].get<0>()==it)
          {
            return j;
          }
      }
    return -1;
  }

/// Set exponent limit for psymbol.
  void symbolic_truncator::set_limit(psym_p it, int n)
  {
    int j=find_expo_limit(it);
    if (j<0)
      {
        std::cout << "Setting new limit." << std::endl;
        vel_.push_back(expo_limit(it,n));
      }
    else
      {
        std::cout << "Modifying existing limit." << std::endl;
        vel_[j].get<1>()=n;
      }
  }
}
