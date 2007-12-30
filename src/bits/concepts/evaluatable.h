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

#ifndef PIRANHA_EVALUATABLE_H
#define PIRANHA_EVALUATABLE_H

#error "Don't use this class, it is just for documentation purposes."

namespace piranha
{
namespace concepts
{
/// Evaluatable class concept.
/**
 * Defines a time-evaluation method. Return value should be decided through the piranha::eval_type type trait.
 * This class serves only for documentation purposes and must not be used.
 */
  template <class Type>
    class evaluatable
  {
    public:
/// Time-evaluation.
/**
 * Evaluate class at time t. Return value can be either double or std::complex<double>, depending on the type
 * of the class.
 * @param[in] t time of evaluation.
 * @see piranha::eval_type for evaluation type trait.
 */
      typename eval_type<Type>::type t_eval(const double &t) const;
  };
}
}
#endif
