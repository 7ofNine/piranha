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

#ifndef PIRANHA_TRUNCATORS_H
#define PIRANHA_TRUNCATORS_H

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../psymbol.h"
#include "../stream_manager.h"

namespace piranha
{
  /// Truncator which does not truncate.
  template <class BaseMultiplier>
    struct no_truncation
  {
    template <class Multiplier>
      no_truncation(const Multiplier &) {}
    template <class Result, class Multiplier>
      bool accept(const Result &, const Multiplier &) const {return true;}
    template <class Cf1, class Cf2, class Key, class Multiplier>
      bool skip(const Cf1 &, const Key &, const Cf2 &, const Key &, const Multiplier &) const {return false;}
  };

   struct __PIRANHA_VISIBLE base_norm_truncator
  {
      static void set(const int &n) throw (unsuitable)
      {
        if (n < 0)
        {
          throw (unsuitable("Please insert a non-negative integer."));
        }
        else if (n == 0)
        {
          m_truncation_level = 0;
        }
        else
        {
          m_truncation_level = std::pow(10.,-n);
        }
      }
      static void print(std::ostream &stream = std::cout);
      static std::string print_to_string();
    private:
      static double m_truncation_level;
  };

  /// Norm-based truncator.
  template <class BaseMultiplier>
    struct norm_truncator:public base_norm_truncator
  {
      template <class Multiplier>
        norm_truncator(const Multiplier &m):m_delta_threshold(
        m.m_s1.calculate_norm(m.m_args_tuple)*m.m_s2.calculate_norm(m.m_args_tuple)*m_truncation_level/
        (2*m.m_s1.template nth_index<0>().size()*m.m_s2.template nth_index<0>().size()))
      {}
      template <class Result, class Multiplier>
        bool accept(const Result &, const Multiplier &) const {return true;}
      template <class Cf1, class Cf2, class Key, class Multiplier>
        bool skip(const Cf1 &c1, const Key &, const Cf2 &c2, const Key &, const Multiplier &m) const
      {
        return (c1.norm(m.m_args_tuple) * c2.norm(m.m_args_tuple) / 2 < m_delta_threshold);
      }
    private:
      const double    m_delta_threshold;
  };

  class __PIRANHA_VISIBLE base_expo_truncator
  {
    protected:
      typedef std::vector<std::pair<psym_p,int> > container_type;
      typedef container_type::iterator iterator;
    public:
      static void limit(const std::string &name, const int &n) throw (not_existing)
      {
        psym_p tmp(psymbol_manager::get_pointer(name));
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
      static void clear(const std::string &name) throw (not_existing)
      {
        psym_p tmp(psymbol_manager::get_pointer(name));
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
        expo_truncator(const Multiplier &m)
      {
        const size_t limits_size = m_expo_limits.size(),
          args_size = m.m_args_tuple.template get<Multiplier::key_type::position>().size();
        for (size_t i = 0; i < limits_size; ++i)
        {
          for (size_t j = 0; j < args_size; ++j)
          {
            if (m_expo_limits[i].first == m.m_args_tuple.template get<Multiplier::key_type::position>()[j])
            {
              m_positions.push_back(std::pair<size_t,int>(j,m_expo_limits[i].second));
              // We can break out, there should not be duplicates inside the arguments list.
              break;
            }
          }
        }
      }
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
      std::vector<std::pair<size_t,int> > m_positions;
      typename BaseMultiplier::key_type   m_tmp_key;
  };
}

#endif
