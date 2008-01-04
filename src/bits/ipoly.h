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

#ifndef PIRANHA_IPOLY_H
#define PIRANHA_IPOLY_H

#include "base_classes/base_ipoly.h"

namespace piranha
{
/// Indexed polynomial.
  template <class Cf, class Index, class Expo>
    class ipoly:public base_ipoly<Cf,Index,Expo,ipoly<Cf,Index,Expo> >
  {
      typedef base_ipoly<Cf,Index,Expo,ipoly<Cf,Index,Expo> > ancestor;
      typedef typename ancestor::im_type im_type;
    public:
      typedef typename ancestor::vector_expo vector_expo;
      ipoly():ancestor::base_ipoly(),private_width_(0)
        {}
      ipoly(const Cf &value):ancestor::base_ipoly(value),private_width_(0)
        {}
      ipoly(const Cf &value, const vector_expo &v):private_width_(v.size())
      {
        ancestor::builder_from_vector(value,v);
      }
      ipoly(const ipoly &p):ancestor::base_ipoly(p),private_width_(p.private_width_)
        {}
      ~ipoly()
        {}
      ipoly &operator=(const ipoly &p)
      {
        if (this != &p)
        {
          private_width_=p.private_width_;
          ancestor::common_assignment(p);
        }
        return *this;
      }
      const uint16 &g_width() const
      {
        return private_width_;
      }
      ipoly &operator+=(const ipoly &p)
      {
        ancestor::addition(p);
        return *this;
      }
      ipoly &operator-=(const ipoly &p)
      {
        ancestor::subtraction(p);
        return *this;
      }
      ipoly &operator*=(const ipoly &p)
      {
        ancestor::mult_by(p);
        return *this;
      }
    private:
      uint16       private_width_;
  };
}

#endif
