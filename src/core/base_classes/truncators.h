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

namespace piranha
{
  /// Truncator which does not truncate.
  template <class Multiplier>
    struct no_truncation
  {
    no_truncation(const Multiplier &) {}
    template <class Term, class ArgsTuple>
      bool accept(const Term &, const ArgsTuple &) const {return true;}
    template <class Cf1, class Cf2, class Key, class ArgsTuple>
      bool skip(const Cf1 &, const Key &, const Cf2 &, const Key &, const ArgsTuple &) const {return false;}
  };

  /// Norm-based truncator.
  template <class Multiplier>
    struct norm_truncator
  {
      norm_truncator(const Multiplier &m):m_delta_threshold(
        m.m_s1.calculate_norm(m.m_args_tuple)*m.m_s2.calculate_norm(m.m_args_tuple)*m.m_s1.get_truncation()/
        (2*m.m_s1.template nth_index<0>().size()*m.m_s2.template nth_index<0>().size())
        )
      {}
      template <class Term, class ArgsTuple>
        bool accept(const Term &, const ArgsTuple &) const {return true;}
      template <class Cf1, class Cf2, class Key, class ArgsTuple>
        bool skip(const Cf1 &c1, const Key &, const Cf2 &c2, const Key &, const ArgsTuple &args_tuple) const
      {
        return (c1.norm(args_tuple) * c2.norm(args_tuple) / 2 < m_delta_threshold);
      }
    private:
      const double    m_delta_threshold;
  };

  /// Truncators for polynomials based on the exponent of one or more variables.
  template <class Series1, class Series2, class ArgsTuple>
    struct poly_exponents_truncator
  {
  };
}

#endif
