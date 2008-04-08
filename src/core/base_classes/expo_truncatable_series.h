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

#ifndef PIRANHA_EXPO_TRUNCATABLE_SERIES_H
#define PIRANHA_EXPO_TRUNCATABLE_SERIES_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../psymbol.h"

namespace piranha
{
  class expo_truncatable_series
  {
      typedef std::vector<std::pair<psym_p,int> > container_type;
      typedef container_type::iterator iterator;
    public:
      static void add_expo_limit(const std::string &name, const int &n)
      {
        std::pair<bool,psym_p> tmp(psymbol_manager::get_pointer(name));
        switch (tmp.first)
        {
          case true:
          {
            iterator it = find_argument(tmp.second);
            if (it == m_expo_limits.end())
            {
              m_expo_limits.push_back(std::pair<psym_p,int>(tmp.second,n));
            }
            else
            {
              it->second = n;
            }
            break;
          }
          case false:
            std::cout << "Failed to add exponent limit, no symbol named \"" << name << "\"." << std::endl;
        }
      }
      static void clear_expo_limits() {m_expo_limits.clear();}
      static void dump_expo_limits()
      {
        const iterator it_f = m_expo_limits.end();
        for (iterator it = m_expo_limits.begin(); it != it_f; ++it)
        {
          std::cout << it->first->name() << ',' << it->second << std::endl;
        }
      }
      static const container_type &get_expo_limits() {return m_expo_limits;}
    private:
      static iterator find_argument(const psym_p &p)
      {
        const iterator it_f = m_expo_limits.end();
        iterator it(m_expo_limits.begin());
        for (; it != it_f; ++it)
        {
          if (it->first == p)
          {
            break;
          }
        }
        return it;
      }
    private:
      static container_type m_expo_limits;
  };
}

#endif
