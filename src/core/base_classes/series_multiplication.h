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

#ifndef PIRANHA_SERIES_MULTIPLICATION_H
#define PIRANHA_SERIES_MULTIPLICATION_H

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_SERIES_MULTIPLICATION_TP_DECL class Derived, template <class, class, class, \
  template <class> class Truncator> class Multiplier, \
  template <class> class Truncator
#define __PIRANHA_SERIES_MULTIPLICATION_TP Derived,Multiplier,Truncator

namespace piranha
{
  template <__PIRANHA_SERIES_MULTIPLICATION_TP_DECL>
    class series_multiplication
  {
    protected:
      // Multiply term-by-term with another series, and place the result into retval.
      // Preconditions:
      // - args_tuple must be the result of a merging of arguments between the two series being multiplied,
      // - retval must be empty.
      template <class Derived2, class ArgsTuple>
        void multiply_by_series(const Derived2 &s2, Derived &retval,
        const ArgsTuple &args_tuple) const
      {
        typedef typename Derived::const_sorted_iterator const_sorted_iterator;
        typedef typename Derived2::const_sorted_iterator const_sorted_iterator2;
        typedef typename Derived::term_type term_type;
        typedef typename Derived2::term_type term_type2;
        // Make sure we are inserting into an empty return value.
        p_assert(retval.template nth_index<0>().empty());
        // Just leave an empty series if this or s2 are zero.
        if (derived_const_cast->template nth_index<0>().empty() or s2.template nth_index<0>().empty())
        {
          ;
        }
        // Optimize if the second series is a pure coefficient series.
        // TODO: test the effectiveness of this by multiplying with single cf series in the first and second place.
        else if (s2.is_single_cf())
        {
          derived_const_cast->multiply_coefficients_by(s2.template nth_index<0>().begin()->m_cf,retval,args_tuple);
        }
        else
        {
          Multiplier<Derived,Derived2,ArgsTuple,Truncator> m(*derived_const_cast,s2,retval,args_tuple);
          m.perform_multiplication();
        }
      }
  };
}

#undef derived_const_cast
#undef derived_cast
#undef __PIRANHA_SERIES_MULTIPLICATION_TP_DECL
#undef __PIRANHA_SERIES_MULTIPLICATION_TP

#endif