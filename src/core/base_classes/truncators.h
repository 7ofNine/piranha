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
  template <class Series1, class Series2, class ArgsTuple>
    class no_truncation
  {
      typedef typename Series1::term_type term_type;
      typedef typename Series1::term_type::cf_type cf_type1;
      typedef typename Series2::term_type::cf_type cf_type2;
      typedef typename Series1::term_type::key_type key_type;
    public:
      no_truncation(const Series1 &, const Series2 &, const ArgsTuple &) {}
      bool accept(const term_type &) const {return true;}
      bool skip(const cf_type1 &, const key_type &, const cf_type2 &, const key_type &) const {return false;}
  };

  /// Norm-based truncator.
  template <class Series1, class Series2, class ArgsTuple>
    class norm_truncator
  {
      typedef typename Series1::term_type term_type;
      typedef typename Series1::term_type::cf_type cf_type1;
      typedef typename Series2::term_type::cf_type cf_type2;
      typedef typename Series1::term_type::key_type key_type;
    public:
      norm_truncator(const Series1 &s1, const Series2 &s2, const ArgsTuple &a):m_args_tuple(a),
        m_delta_threshold(
        s1.calculate_norm(m_args_tuple)*s2.calculate_norm(m_args_tuple)*s1.get_truncation()/
        (2*s1.template nth_index<0>().size()*s2.template nth_index<0>().size())
        )
      {}
      bool accept(const term_type &) const {return true;}
      bool skip(const cf_type1 &c1, const key_type &, const cf_type2 &c2, const key_type &) const
      {
        return (c1.norm(m_args_tuple) * c2.norm(m_args_tuple) / 2 < m_delta_threshold);
      }
    private:
      const ArgsTuple &m_args_tuple;
      const double    m_delta_threshold;
  };

  /// Truncators for polynomials based on the exponent of one or more variables.
  template <class Series1, class Series2, class ArgsTuple>
    class poly_exponents_truncator
  {
      typedef typename Series1::term_type term_type;
      typedef typename Series1::term_type::cf_type cf_type1;
      typedef typename Series2::term_type::cf_type cf_type2;
      typedef typename Series1::term_type::key_type key_type;
    public:
  };
}

#endif
