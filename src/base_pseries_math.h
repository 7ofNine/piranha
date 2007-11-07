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
/// Basic assignment.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::basic_assignment(const
    base_pseries &ps2)
  {
    if (this==&ps2)
    {
      return;
    }
    *s_series_set()=*ps2.g_series_set();
    cf_s_vec_=ps2.cf_s_vec_;
    trig_s_vec_=ps2.trig_s_vec_;
    lin_args_=ps2.lin_args_;
    static_cast<Derived *>(this)->assignment_hook(ps2);
    std::cout << "Assignment operator!" << std::endl;
  }

/// Assignment from series with differen coefficient.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::generic_series_assignment(const Derived2 &ps2)
  {
    if ((void *)this==(void *)&ps2)
    {
      return;
    }
    typedef typename Derived2::it_s_index it_s_index2;
    s_series_set()->clear();
    cf_s_vec_=ps2.cf_s_vec();
    trig_s_vec_=ps2.trig_s_vec();
    lin_args_=ps2.lin_args();
    const it_s_index2 it_f=ps2.g_s_index().end();
// TODO: use hinted insertion.
    for (it_s_index2 it=ps2.g_s_index().begin();it!=it_f;++it)
    {
      insert(*it);
    }
    static_cast<Derived *>(this)->assignment_hook(ps2);
    std::cout << "Generic assignment operator!" << std::endl;
  }

// Base merge operator
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::alg_sum_lin_args(const Derived2 &ps2,
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

  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::merge_with(const Derived2 &ps2, bool sign)
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
    it_s_index it_hint=g_s_index().end();
// NOTE: At this point this' size is greater or equal to ps2'
    for (typename Derived2::ancestor::it_h_index it=ps2.g_h_index().begin();
      it!=ps2.g_h_index().end();++it)
    {
      it_hint=insert(*it,sign,&it_hint);
    }
  }

// Merge with a generic entity - NOT with another series
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::generic_merge(const T &x)
  {
// Build a series from x
    Derived tmp=Derived(cf_type(x),*static_cast<Derived *>(this));
// Merge with this
    merge_with(tmp);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Cf2, class LightTermPair>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::term_by_term_multiplication_trig(
    const term_type &t1, const Term<Cf2,trig_type> &t2, LightTermPair &term_pair,
    cf_type &new_c)
  {
    t1.g_trig()->trigmult(*t2.g_trig(),*term_pair.template get
      <0>().s_trig(),
      *term_pair.template get<1>().s_trig());
    *term_pair.template get
      <0>().s_cf()=*term_pair.template get
      <1>().s_cf()=new_c;
    if (t1.g_flavour()==t2.g_flavour())
    {
      term_pair.template get
        <0>().s_trig()->s_flavour()=term_pair.template get
        <1>().s_trig()->s_flavour()=true;
      if(!t1.g_flavour())
      {
        *term_pair.template get
          <1>().s_cf()*=-1;
      }
    }
    else
    {
      term_pair.template get
        <0>().s_trig()->s_flavour()=term_pair.template get
        <1>().s_trig()->s_flavour()=false;
      if(t1.g_flavour())
      {
        *term_pair.template get
            <0>().s_cf()*=-1;
      }
    }
  }

/// Series multiplication preliminaries.
/**
 * Perform some preliminary activity regarding series multiplication. If some checks are successful and arguments
 * were merged successfully return true, otherwise return false.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline bool base_pseries<__PIRANHA_BASE_PS_TP>::series_multiplication_preliminaries
    (const Derived2 &ps2, Derived &retval)
  {
    if (length()==0 || ps2.length()==0)
    {
      std::cout << "Zero stuff." << std::endl;
      return false;
    }
    if (!math::is_zero_vec(lin_args())||!math::is_zero_vec(ps2.lin_args()))
    {
      std::cout << "Non-zero linargs in multiplication." << std::endl;
      std::exit(1);
      return false;
    }
    if (!merge_args(ps2))
    {
      std::cout << "Args are not compatible during multiplication, returning null series." << std::endl;
      std::exit(1);
      return false;
    }
    action_assert(retval.merge_args(*static_cast<Derived *>(this)));
    return true;
  }

/// Optimize series multiplication for simple series.
/**
 * If at least one of two multiplied series is formed by a single term with trigonometric part null and cosine,
 * perform a cheaper coefficient multiplication instead of a term multiplication.
 * Returns true if such optimization could be performed, false otherwise.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline bool base_pseries<__PIRANHA_BASE_PS_TP>::series_multiplication_optimize_for_cf
    (const Derived2 &ps2)
  {
    if (ps2.is_cf())
    {
      std::cout << "Cf1\n";
      static_cast<Derived *>(this)->cf_multiplication(*ps2.g_s_index().begin()->g_cf());
      return true;
    }
    else if (is_cf())
    {
      cf_type tmp(*g_s_index().begin()->g_cf());
      generic_series_assignment(ps2);
      static_cast<Derived *>(this)->cf_multiplication(tmp);
      std::cout << "Cf2\n";
      return true;
    }
    return false;
  }

/// Series multiplication.
/**
 * Requires some methods to be implemented in derived classes.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::series_multiplication(const Derived2 &ps2)
  {
    Derived *derived_cast=static_cast<Derived *>(this);
    Derived retval;
    if (series_multiplication_preliminaries(ps2,retval))
    {
      if (series_multiplication_optimize_for_cf(ps2))
      {
        return;
      }
      derived_cast->multiply_terms(ps2,retval);
    }
    swap(retval);
  }

/// Multiplication by a generic entity.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::generic_multiplication(const T &c)
  {
    if (empty())
    {
      return;
    }
    if (!math::is_zero_vec(lin_args()))
    {
      std::cout << "Non-zero linargs in *= T!" << std::endl;
      std::exit(1);
    }
    Derived tmp_ps;
    tmp_ps.merge_args(*static_cast<Derived *>(this));
    term_type tmp_term;
    it_s_index it_hint=tmp_ps.g_s_index().end();
    const it_s_index it_f=g_s_index().end();
    for (it_s_index it=g_s_index().begin();it!=it_f;++it)
    {
      tmp_term=*it;
      *tmp_term.s_cf()*=c;
      it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
    }
    swap(tmp_ps);
  }

/// Generic division.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::generic_division(const T &x)
  {
    if (x==0)
    {
      std::cout << "ERROR: division by zero in /=, returning self." << std::endl;
      std::abort();
      return;
    }
    if (empty())
    {
      return;
    }
    if (!math::is_zero_vec(lin_args()))
    {
// NOTICE: maybe here we could deal with exact int/int divisions. Not really important
// ATM though.
      std::cout << "Non-zero linargs in /= int!" << std::endl;
      std::exit(1);
    }
    Derived tmp_ps;
    tmp_ps.merge_args(*static_cast<Derived *>(this));
    term_type tmp_term;
    it_s_index it_hint=tmp_ps.g_s_index().end();
    const it_s_index it_f=g_s_index().end();
    for (it_s_index it=g_s_index().begin();it!=it_f;++it)
    {
      tmp_term=*it;
      *tmp_term.s_cf()/=x;
      it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
    }
    swap(tmp_ps);
  }

/// Multiplication by an integer.
/**
 * This is specialized because we have to take care of lin_args.
 */
/// Generic division.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::mult_by_int(int n)
  {
    const vector_int16 old_lin_args=lin_args();
    size_t j;
    const size_t w=lin_args().size();
// Zero the linargs, otherwise the generic *= operator complains
    for (j=0;j<w;++j)
    {
      lin_args()[j]=0;
    }
// Now perform the generic multiplication
    generic_multiplication(n);
// Multiply the old linargs and restore them
    for (j=0;j<w;++j)
    {
      lin_args()[j]=old_lin_args[j]*n;
    }
  }
}
#endif
