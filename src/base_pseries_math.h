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
    base_pseries<Cf, Trig, Term, I, Derived> &ps2)
  {
    if (this==&ps2)
    {
      return;
    }
    set_=ps2.set_;
    norm_=ps2.norm_;
    cf_s_vec_=ps2.cf_s_vec_;
    trig_s_vec_=ps2.trig_s_vec_;
    lin_args_=ps2.lin_args_;
    std::cout << "Assignment operator!" << std::endl;
  }

/************************/
/* Low level operations */
/************************/

// Base merge operator
// -------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2, class Derived2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::alg_sum_lin_args(const base_pseries<Cf2, trig_type, Term, I, Derived2> &ps2,
    bool sign)
  {
    vector_mult_t tmp(trig_width());
    if (sign)
    {
      vec_add(lin_args_,ps2.lin_args(),ps2.lin_args(),lin_args_,tmp);
    }
    else
    {
      vec_sub(lin_args_,ps2.lin_args(),ps2.lin_args(),lin_args_,tmp);
    }
    lin_args_=tmp;
  }

  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2, class Derived2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::merge_with(const base_pseries<Cf2, trig_type, Term, I, Derived2> &ps2, bool sign)
  {
    if ((void *)&ps2==(void *)this)
    {
      if (sign)
      {
        base_pseries tmp_ps(*this);
        tmp_ps.merge_with(ps2,sign);
        swap(tmp_ps);
      }
      else
      {
        base_pseries tmp_ps;
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
    for (typename base_pseries<Cf2, trig_type, Term, I, Derived2>::it_h_index it=ps2.h_index().begin();
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
    base_pseries tmp=base_pseries(cf_type(x));
// Merge with this
    merge_with(tmp);
  }


// Multiplication by a generic entity
// ----------------------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class T>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::generic_multiplication(const T &c)
  {
    if (length()==0)
    {
      return;
    }
    if (!is_zero_vec(lin_args_))
    {
      std::cout << "Non-zero linargs in *= T!" << std::endl;
      std::exit(1);
    }
    base_pseries tmp_ps;
    tmp_ps.merge_args(*this);
    term_type tmp_term;
    it_s_index it_hint=tmp_ps.s_index().end();
    for (it_s_index it=s_index().begin();it!=s_index().end();++it)
    {
      tmp_term=*it;
      tmp_term.s_cf()*=c;
      it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
    }
    swap(tmp_ps);
  }

/// Division by an integer.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Integer>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::basic_div_by_int(const Integer &n)
  {
    if (n==0)
    {
      std::cout << "ERROR: division by zero in /= int, returning self." << std::endl;
      std::abort();
      return;
    }
    if (length()==0)
    {
      return;
    }
    if (!is_zero_vec(lin_args_))
    {
// NOTICE: maybe here we could deal with exact int/int divisions. Not really important
// ATM though.
      std::cout << "Non-zero linargs in /= int!" << std::endl;
      std::exit(1);
    }
    base_pseries tmp_ps;
    tmp_ps.merge_args(*this);
    term_type tmp_term;
    it_s_index it_hint=tmp_ps.s_index().end();
    for (it_s_index it=s_index().begin();it!=s_index().end();++it)
    {
      tmp_term=*it;
      tmp_term.s_cf()/=n;
      it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
    }
    swap(tmp_ps);
  }

// Multiplication by an integer
// ----------------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::mult_by_int(int n)
  {
    const vector_mult_t old_lin_args=lin_args_;
    unsigned int j;
// Zero the linargs, otherwise the generic *= operator complains
    for (j=0;j<lin_args_.size();++j)
    {
      lin_args_[j]=0;
    }
// Now perform the generic multiplication
    generic_multiplication(n);
// Multiply the old linargs and restore them
    const unsigned int w=lin_args_.size();
    for (j=0;j<w;++j)
    {
      lin_args_[j]=old_lin_args[j]*n;
    }
  }
}
#endif
