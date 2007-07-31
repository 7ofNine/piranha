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

#ifndef PIRANHA_BASE_PSERIES_MATH_H
#define PIRANHA_BASE_PSERIES_MATH_H

namespace piranha
{
// Assignment operator
// -------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::basic_assignment(const
    base_pseries &ps2)
  {
    if (this==&ps2)
    {
      return;
    }
    set_=ps2.set_;
    cf_s_vec_=ps2.cf_s_vec_;
    trig_s_vec_=ps2.trig_s_vec_;
    lin_args_=ps2.lin_args_;
    static_cast<Derived *>(this)->assignment_hook(ps2);
    std::cout << "Assignment operator!" << std::endl;
  }

/************************/
/* Low level operations */
/************************/

// Base merge operator
// -------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Derived2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::alg_sum_lin_args(const Derived2 &ps2,
    bool sign)
  {
    if (sign)
    {
      math::vec_add(lin_args_,ps2.lin_args(),lin_args_);
    }
    else
    {
      math::vec_sub(lin_args_,ps2.lin_args(),lin_args_);
    }
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Derived2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::merge_with(const Derived2 &ps2, bool sign)
  {
    if ((void *)&ps2==(void *)this)
    {
      if (sign)
      {
        Derived tmp_ps(*static_cast<Derived *>(this));
        tmp_ps.merge_with(ps2,sign);
        swap(tmp_ps);
      }
      else
      {
        Derived tmp_ps;
        tmp_ps.merge_args(*this);
        tmp_ps.lin_args_=lin_args();
        swap(tmp_ps);
      }
      return;
    }
// Check that trig_args are compatible
    if (!merge_args(ps2))
    {
      std::cout << "trig_args are not compatible, returning self." << std::endl;
      std::exit(1);
      return;
    }
// Sum/sub lin_args
    alg_sum_lin_args(ps2,sign);
// Use hint, since as we add terms we have an idea of where they are going to be placed
    it_s_index it_hint=s_index().end();
// NOTE: At this point this' size is greater or equal to ps2'
    for (typename Derived2::ancestor::it_h_index it=ps2.h_index().begin();
      it!=ps2.h_index().end();++it)
    {
      it_hint=insert(*it,sign,&it_hint);
    }
  }

// Merge with a generic entity - NOT with another series
// -----------------------------------------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class T>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::generic_merge(const T &x)
  {
// Build a series from x
    Derived tmp=Derived(cf_type(x),*static_cast<Derived *>(this));
// Merge with this
    merge_with(tmp);
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::term_by_term_multiplication_trig(
    const term_type &t1, const Term<Cf2,trig_type> &t2, boost::tuple<term_type &,term_type &> &term_pair,
    cf_type &new_c)
  {
    if (t1.g_flavour())
    {
      if(t2.g_flavour())
      {
        t1.g_trig().trigmult(t2.g_trig(),term_pair.template get
          <0>().s_trig(),
          term_pair.template get<1>().s_trig());
        term_pair.template get
          <0>().s_cf()=term_pair.template get
          <1>().s_cf()=new_c;
        term_pair.template get
          <0>().s_flavour()=term_pair.template get
          <1>().s_flavour()=true;
      }
      else
      {
        t1.g_trig().trigmult(t2.g_trig(),term_pair.template get
          <0>().s_trig(),
          term_pair.template get<1>().s_trig());
        term_pair.template get
          <0>().s_cf()=-new_c;
        term_pair.template get
          <1>().s_cf()=new_c;
        term_pair.template get
          <0>().s_flavour()=term_pair.template get
          <1>().s_flavour()=false;
      }
    }
    else
    {
      if(t2.g_flavour())
      {
        t1.g_trig().trigmult(t2.g_trig(),term_pair.template get
          <0>().s_trig(),
          term_pair.template get<1>().s_trig());
        term_pair.template get
          <0>().s_cf()=term_pair.template get
          <1>().s_cf()=new_c;
        term_pair.template get
          <0>().s_flavour()=term_pair.template get
          <1>().s_flavour()=false;
      }
      else
      {
        t1.g_trig().trigmult(t2.g_trig(),term_pair.template get
          <0>().s_trig(),
          term_pair.template get<1>().s_trig());
        term_pair.template get
          <0>().s_cf()=new_c;
        term_pair.template get
          <1>().s_cf()=-new_c;
        term_pair.template get
          <0>().s_flavour()=term_pair.template get
          <1>().s_flavour()=true;
      }
    }
  }
}

#endif
