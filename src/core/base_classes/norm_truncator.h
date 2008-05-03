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

#ifndef PIRANHA_NORM_TRUNCATOR_H
#define PIRANHA_NORM_TRUNCATOR_H

#include <cmath>
#include <iostream>
#include <string>

#include "../config.h"
#include "../exceptions.h"

namespace piranha
{
  struct __PIRANHA_VISIBLE base_norm_truncator
  {
      static void set(const int &n)
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
        m.m_s1.b_norm(m.m_args_tuple)*m.m_s2.b_norm(m.m_args_tuple)*m_truncation_level/
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
      const double  m_delta_threshold;
  };
}

#endif
