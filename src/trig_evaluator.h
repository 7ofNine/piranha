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

#ifndef PIRANHA_TRIG_EVALUATOR_H
#define PIRANHA_TRIG_EVALUATOR_H

#include <boost/tuple/tuple.hpp>
#include <vector>

#include "common_typedefs.h"
#include "p_assert.h"

namespace piranha
{
  template <class DerivedPs>
    class trig_evaluator
  {
      typedef std::vector<complex_double> vector_complex;
      typedef boost::tuple<vector_complex,vector_complex> internal_tuple;
      typedef std::vector<internal_tuple> container_type;
    public:
      trig_evaluator(const DerivedPs *ps, const double &t):private_ps_(ps),private_width_(ps->trig_width()),
        private_value_(t),private_container_(private_width_)
      {
// Initialize powers 1 and -1 for all trigonometric arguments.
        for (size_t i=0;i<private_width_;++i)
        {
          private_container_[i].get<0>().push_back(std::polar(1.,private_ps_->trig_s_vec()[i]->t_eval(t)));
          private_container_[i].get<1>().push_back(complex_double(1.)/
            private_container_[i].get<0>()[0]);
        }
      }
      complex_double request_complexp(const size_t &index, const int &power_)
      {
        p_assert(power_!=0);
// Make sure we are not going outside container's boundaries.
        p_assert(index < private_width_);
        int power=power_;
        vector_complex *exp_vec=&private_container_[index].get<0>();
// Change sign and pointer to vector of complex exponentials if negative.
        if (power_<0)
        {
          power=-power_;
          exp_vec=&private_container_[index].get<1>();
        }
// Add the missing elements.
        while (exp_vec->size() < (size_t)power)
        {
          exp_vec->push_back((*exp_vec)[0]*((*exp_vec)[exp_vec->size()-1]));
        }
        return (*exp_vec)[(size_t)power-1];
      }
      const double &value() const
      {
        return private_value_;
      }
      const DerivedPs *ps() const
      {
        return private_ps_;
      }
      const size_t &width() const
      {
        return private_width_;
      }
    private:
// Data members.
      DerivedPs const   *private_ps_;
      const size_t      private_width_;
      const double      private_value_;
      container_type    private_container_;
  };
}

#endif
