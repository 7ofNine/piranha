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

#ifndef PIRANHA_TYPE_TRAITS_H
#define PIRANHA_TYPE_TRAITS_H

#include <complex>

namespace piranha
{
  /// Evaluation type trait.
  /**
   * Specifies the type of evaluation: default is double.
   */
  template <class T>
    struct eval_type
  {
    typedef double type;
  };

  /// Complex specialization for evaluation type trait.
  /**
   * Evaluation type is the complex counterpart of real evaluation type.
   */
  template <class T>
    struct eval_type<std::complex<T> >
  {
    typedef std::complex<typename eval_type<T>::type> type;
  };

  /// Representation used for coefficients and keys of piranha::base_series objects during multiplication.
  template <class T>
    struct series_mult_rep
  {
    /// Representation type.
    typedef T type;
    /// Get a const reference to the element represented by x.
    static const T &get(const type &x) {return x;}
    /// Assign T value source to the representation type res.
    static void assign(type &res, const T &source) {res=source;}
  };
}
#endif
