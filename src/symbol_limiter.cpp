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

#include "symbol_limiter.h"

namespace piranha
{
/// Initialization of static member of piranha::symbol_limiter.
  symbol_limiter::limits_map symbol_limiter::lmap_;

/// Set exponent limit for psymbol.
  void symbol_limiter::set_limit(psym_p it, int n)
  {
    if (it==psymbol_manager::end())
      {
        std::cout << "Won't set limit for invalid psymbol pointer." << std::endl;
        std::abort();
        return;
      }
    if (find_expo_limit(it)==lmap_.end())
      {
        std::cout << "Setting new limit." << std::endl;
        lmap_[it]=n;
      }
    else
      {
        std::cout << "Modifying existing limit." << std::endl;
        mit->second=n;
      }
  }

/// Clear limits.
  void symbol_limiter::clear()
  {
    lmap_.clear();
  }

/// Print to screen limits list.
  void symbol_limiter::put()
  {
    for (map_iterator it=lmap_.begin();it!=lmap_.end();++it)
      {
        std::cout << "Limit for '" << it->first->name() << "': " << it->second << std::endl;
      }
  }
}
