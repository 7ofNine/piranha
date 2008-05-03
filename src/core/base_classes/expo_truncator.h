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

#ifndef PIRANHA_EXPO_TRUNCATOR_H
#define PIRANHA_EXPO_TRUNCATOR_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../exceptions.h"
#include "../psym.h"

namespace piranha
{
  /// Base exponent truncator.
  /**
   * Internally it stores a static list of piranha::psym_p + int pairs representing the limits
   * for the exponents of symbols. This class can be used by those classes that need the
   * exponent limits (such as truncators for polynomial multiplication, classes that perform power series
   * expansions, etc.).
   */
  class __PIRANHA_VISIBLE base_expo_truncator
  {
    protected:
      typedef std::vector<std::pair<psym_p,int> > container_type;
      typedef container_type::iterator iterator;
    public:
      static void limit(const std::string &name, const int &n)
      {
        psym_p tmp(psym_manager::get_pointer(name));
        iterator it = find_argument(tmp);
        if (it == m_expo_limits.end())
        {
          m_expo_limits.push_back(std::pair<psym_p,int>(tmp,n));
        }
        else
        {
          it->second = n;
        }
      }
      static void clear_all();
      static void print(std::ostream &stream = std::cout);
      static std::string print_to_string();
      static void clear(const std::string &name)
      {
        psym_p tmp(psym_manager::get_pointer(name));
        iterator it = find_argument(tmp);
        if (it == m_expo_limits.end())
        {
          throw (not_existing(std::string("Symbol ") +"\"" + name + "\" does not have an exponent limit set."));
        }
        else
        {
          m_expo_limits.erase(it);
        }
      }
      /// Transform the list of psymbol limits into a list of positions - limits given a piranha::vector_psym_p.
      static std::vector<std::pair<size_t,int> > get_positions_limits(const vector_psym_p &v)
      {
        std::vector<std::pair<size_t,int> > retval;
        const size_t limits_size = m_expo_limits.size(),
          args_size = v.size();
        for (size_t i = 0; i < limits_size; ++i)
        {
          for (size_t j = 0; j < args_size; ++j)
          {
            if (m_expo_limits[i].first == v[j])
            {
              retval.push_back(std::pair<size_t,int>(j,m_expo_limits[i].second));
              // We can break out, there should not be duplicates inside the arguments list.
              break;
            }
          }
        }
        return retval;
      }
    private:
      static iterator find_argument(const psym_p &);
    protected:
      static container_type m_expo_limits;
  };

  /// Truncators for polynomials based on the exponent of one or more variables.
  template <class BaseMultiplier>
    class expo_truncator:public base_expo_truncator
  {
    public:
      template <class Multiplier>
        expo_truncator(const Multiplier &m):
        m_positions(base_expo_truncator::get_positions_limits(m.m_args_tuple.template get<Multiplier::key_type::position>()))
      {}
      template <class Multiplier>
        bool accept(const max_fast_int &n, const Multiplier &m)
      {
        switch (m_positions.size() == 0)
        {
          case true:
            return true;
          default:
            m.decode(m_tmp_key,n);
            return m_tmp_key.test_expo_limits(m_positions);
        }
      }
      template <class Term, class Multiplier>
        bool accept(const Term &t, const Multiplier &) const
      {
        switch (m_positions.size() == 0)
        {
          case true:
            return true;
          default:
            return t.m_key.test_expo_limits(m_positions);
        }
      }
      template <class Cf1, class Cf2, class Key, class Multiplier>
        bool skip(const Cf1 &, const Key &, const Cf2 &, const Key &, const Multiplier &) const
      {
        return false;
      }
    private:
      const std::vector<std::pair<size_t,int> > m_positions;
      typename BaseMultiplier::key_type         m_tmp_key;
  };
}

#endif
