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

#ifndef PIRANHA_PS_TERM_MATH_H
#define PIRANHA_PS_TERM_MATH_H

namespace piranha
{
/// Assignment operator.
  template <class Cf,class Trig>
    inline ps_term<Cf,Trig> &ps_term<Cf,Trig>::operator=(const ps_term<Cf,Trig> &t2)
  {
    if (this==&t2)
    {
      return *this;
    }
    c_=t2.c_;
    trig_args_=t2.trig_args_;
    flavour_=t2.flavour_;
    return *this;
  }

  template <class Cf,class Trig>
    template <class T,class U>
    inline void ps_term<Cf,Trig>::mult_by(const T &t2, boost::tuple<U,U> &term_pair) const
  {
    cf_type new_c=c_;
    new_c*=t2.c();
    new_c/=2;
    if (flavour_)
    {
      if(t2.flavour())
      {
        trig_args_.trigmult(t2.trig_args(),term_pair.template get
          <0>().trig_args(),
          term_pair.template get<1>().trig_args());
        term_pair.template get
          <0>().c()=term_pair.template get
          <1>().c()=new_c;
        term_pair.template get
          <0>().flavour()=term_pair.template get
          <1>().flavour()=true;
      }
      else
      {
        trig_args_.trigmult(t2.trig_args(),term_pair.template get
          <0>().trig_args(),
          term_pair.template get<1>().trig_args());
        term_pair.template get
          <0>().c()=-new_c;
        term_pair.template get
          <1>().c()=new_c;
        term_pair.template get
          <0>().flavour()=term_pair.template get
          <1>().flavour()=false;
      }
    }
    else
    {
      if(t2.flavour())
      {
        trig_args_.trigmult(t2.trig_args(),term_pair.template get
          <0>().trig_args(),
          term_pair.template get<1>().trig_args());
        term_pair.template get
          <0>().c()=term_pair.template get
          <1>().c()=new_c;
        term_pair.template get
          <0>().flavour_=term_pair.template get
          <1>().flavour()=false;
      }
      else
      {
        trig_args_.trigmult(t2.trig_args(),term_pair.template get
          <0>().trig_args(),
          term_pair.template get<1>().trig_args());
        term_pair.template get
          <0>().c()=new_c;
        term_pair.template get
          <1>().c()=-new_c;
        term_pair.template get
          <0>().flavour()=term_pair.template get
          <1>().flavour()=true;
      }
    }
  }
}
#endif
