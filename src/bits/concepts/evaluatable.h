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

#include "boost/static_assert.hpp"

namespace piranha
{
namespace concepts
{
/// Evaluatable class concept.
/**
 * Defines a typedef and a time-evaluation method. This class serves only for documentation purposes.
 */
  template <Class EvalType>
    class evaluatable
  {
    public:
/// Evaluation type.
/**
 * This will be the result of the evaluation method.
 */
      typedef EvalType eval_type;
/// Time-evaluation.
/**
 * Evaluate class at time t.
 * @param[in] t time of evaluation.
 */
      eval_type t_eval(const double &t) const
      {
        BOOST_STATIC_ASSERT(sizeof(EvalType) == 0);
        (void)t;
        return eval_type();
      }
    private:
      evaluatable() {BOOST_STATIC_ASSERT(sizeof(EvalType) == 0);}
  };
}
}
#endif
