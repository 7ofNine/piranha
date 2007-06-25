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

// Low-level mutiplication of terms
// --------------------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2,class Derived2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::mult_terms(const base_pseries<Cf2, trig_type, Term, I, Derived2> &ps2,
    base_pseries &retval, const double &Delta)
  {
    const double Delta_threshold=Delta/(2*length()*ps2.length());
    size_t n=0;
// NOTE: at this point retval's width() is greater or equal to _both_ this
// and ps2. It's the max width indeed.
    p_assert(math::max(trig_width(),ps2.trig_width())==retval.trig_width());
    term_type tmp1, tmp2;
    boost::tuple<term_type &, term_type &> term_pair(tmp1,tmp2);
    const it_s_index it1_f=s_index().end();
    const typename base_pseries<Cf2, trig_type, Term, I, Derived2>::it_s_index it2_f=ps2.s_index().end();
    typename base_pseries<Cf2, trig_type, Term, I, Derived2>::it_s_index it2;
    it_s_index it1, it_hint=retval.s_index().end();
    for (it1=s_index().begin();it1!=it1_f;++it1)
    {
      it2=ps2.s_index().begin();
      if ((it1->norm(cf_s_vec_)*it2->norm(ps2.cf_s_vec()))/2<Delta_threshold)
      {
        break;
      }
      for (;it2!=it2_f;++it2)
      {
// We are going to calculate a term's norm twice... We need to profile
// this at a later stage and see if it is worth to store the norm inside
// the term.
        if ((it1->norm(cf_s_vec_)*it2->norm(ps2.cf_s_vec()))/2<Delta_threshold)
        {
          break;
        }
        it1->mult_by(*it2,term_pair);
// Before insertion we change the sign of trigonometric parts if necessary.
// This way we won't do a copy inside insertion function.
        if (term_pair.template get
          <0>().g_trig().sign()<0)
        {
          term_pair.template get<0>().invert_trig_args();
        }
        if (term_pair.template get
          <1>().g_trig().sign()<0)
        {
          term_pair.template get<1>().invert_trig_args();
        }
        it_hint=retval.insert(term_pair.template get
          <0>(),true,&it_hint);
        it_hint=retval.insert(term_pair.template get
          <1>(),true,&it_hint);
        ++n;
      }
    }
//retval.cumulative_crop(Delta);
    std::cout << "w/o trunc=" << length()*ps2.length() << "\tw/ trunc=" << n << std::endl;
    std::cout << "Out length=" << retval.length() << std::endl;
  }

// Basic multiplication
// --------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class Cf2, class Derived2>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::basic_ps_mult(const
    base_pseries<Cf2, Trig, Term, I, Derived2> &ps2)
  {
    base_pseries retval;
// OPTIMIZE: optimize in the case one series is a c value
// If one length is zero do not do anything
    if (length()!=0 && ps2.length()!=0)
    {
      if (!is_zero_vec(lin_args_)||!is_zero_vec(ps2.lin_args()))
      {
        std::cout << "Non-zero linargs!" << std::endl;
        std::exit(1);
      }
      if (!merge_args(ps2))
      {
        std::cout << "args are not compatible, returning self." << std::endl;
        std::exit(1);
        return;
      }
      p_assert(retval.merge_args(*this));
      const double Delta=norm()*ps2.norm()*settings_manager::prec();
      mult_terms(ps2,retval,Delta);
    }
    else
    {
      std::cout << "Zero stuff" << std::endl;
    }
    swap(retval);
  }

// Multiplication by a generic entity
// ----------------------------------
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    template <class T>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::generic_mult(const T &c)
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
      tmp_term.s_c()*=c;
      it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
    }
    swap(tmp_ps);
  }

/// Division by an integer.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::basic_div_by_int(int n)
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
      tmp_term.s_c()/=n;
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
    generic_mult(n);
// Multiply the old linargs and restore them
    const unsigned int w=lin_args_.size();
    for (j=0;j<w;++j)
    {
      lin_args_[j]=old_lin_args[j]*n;
    }
  }
}
#endif
