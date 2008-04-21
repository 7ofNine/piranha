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

#ifndef PIRANHA_TRIG_EVALUATOR_H
#define PIRANHA_TRIG_EVALUATOR_H

#include <complex>
#include <utility>
#include <vector>

#include "../p_assert.h"

namespace piranha
{
  template <class DerivedPs>
    class trig_evaluator
  {
      typedef std::vector<std::complex<double> > vector_complex;
      typedef std::pair<vector_complex,vector_complex> internal_pair;
      typedef std::vector<internal_pair> container_type;
    public:
      trig_evaluator(const DerivedPs *ps, const double &t):private_ps_(ps),private_width_(ps->trig_width()),
        private_value_(t),private_container_(private_width_)
      {
        // Initialize powers 1 and -1 for all trigonometric arguments.
        for (size_t i=0;i<private_width_;++i)
        {
          private_container_[i].first.push_back(std::polar(1.,private_ps_->arguments().template get<1>()[i]->t_eval(t)));
          private_container_[i].second.push_back(std::complex<double>(1.)/
            private_container_[i].first[0]);
        }
      }
      std::complex<double> request_complexp(const size_t &index, const int &power_)
      {
        p_assert(power_!=0);
        // Make sure we are not going outside container's boundaries.
        p_assert(index < private_width_);
        int power=power_;
        vector_complex *exp_vec=&private_container_[index].first;
        // Change sign and pointer to vector of complex exponentials if negative.
        if (power_<0)
        {
          power=-power_;
          exp_vec=&private_container_[index].second;
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
