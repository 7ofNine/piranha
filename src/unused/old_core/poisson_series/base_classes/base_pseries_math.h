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

#ifndef PIRANHA_BASE_PSERIES_MATH_H
#define PIRANHA_BASE_PSERIES_MATH_H

#include "base_pseries_ta_macros.h"

namespace piranha
{
  /// Assign self.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::assign(const
    Derived &ps2)
  {
    if (static_cast<Derived *>(this) == &ps2)
    {
      return *static_cast<Derived *>(this);
    }
    *s_series_set()=*ps2.g_series_set();
    arguments().template get<0>()=ps2.arguments().template get<0>();
    arguments().template get<1>()=ps2.arguments().template get<1>();
    lin_args()=ps2.lin_args();
    static_cast<Derived *>(this)->assignment_hook(ps2);
    std::cout << "Assignment operator!" << std::endl;
    return *static_cast<Derived *>(this);
  }

  /// Assignment from different series.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::assign_series(const Derived2 &ps2)
  {
    // TODO: use static_cast.
    if ((void *)this==(void *)&ps2)
    {
      return *static_cast<Derived *>(this);
    }
    typedef typename Derived2::it_s_index it_s_index2;
    s_series_set()->clear();
    arguments().template get<0>()=ps2.arguments().template get<0>();
    arguments().template get<1>()=ps2.arguments().template get<1>();
    lin_args()=ps2.lin_args();
    const it_s_index2 it_f=ps2.g_s_index().end();
    it_s_index it_hint = g_s_index().end();
    for (it_s_index2 it=ps2.g_s_index().begin();it!=it_f;++it)
    {
      it_hint=insert(*it,it_hint);
    }
    static_cast<Derived *>(this)->assignment_hook(ps2);
    std::cout << "Generic assignment operator!" << std::endl;
    return *static_cast<Derived *>(this);
  }

  /// Assignment from generic type.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::assign_generic(const T &x)
  {
    s_series_set()->clear();
    arguments().template get<0>().clear();
    arguments().template get<1>().clear();
    lin_args().clear();
    insert(term_type(cf_type(x),trig_type()));
    std::cout << "Generic assignment operator!" << std::endl;
    return *static_cast<Derived *>(this);
  }

  /// Assign int.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::assign(int n)
  {
    return assign_generic(n);
  }

  /// Assign double.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::assign(const double &x)
  {
    return assign_generic(x);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2, bool Sign>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::alg_sum_lin_args(const Derived2 &ps2)
  {
    if (Sign)
    {
      math::vec_add(lin_args(),ps2.lin_args(),lin_args());
    }
    else
    {
      math::vec_sub(lin_args(),ps2.lin_args(),lin_args());
    }
  }

  // Base merge operator
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2, bool Sign>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::merge_with_series(const Derived2 &ps2)
  {
    if ((void *)&ps2==(void *)this)
    {
      if (Sign)
      {
        Derived tmp_ps(*static_cast<Derived *>(this));
        tmp_ps.merge_with_series<Derived2,Sign>(ps2);
        swap(tmp_ps);
      }
      else
      {
        Derived tmp_ps;
        tmp_ps.merge_args(*this);
        tmp_ps.lin_args()=lin_args();
        swap(tmp_ps);
      }
      return *static_cast<Derived *>(this);
    }
    try
    {
      merge_args(ps2);
      // Sum/sub lin_args
      alg_sum_lin_args<Derived2,Sign>(ps2);
      // Use hint, since as we add terms we have an idea of where they are going to be placed
      it_s_index it_hint=g_s_index().end();
      // NOTE: At this point this' size is greater or equal to ps2'
      for (typename Derived2::ancestor::it_h_index it=ps2.g_h_index().begin();
        it!=ps2.g_h_index().end();++it)
      {
        it_hint=ancestor::template insert<true,Sign>(*it,it_hint);
      }
      return *static_cast<Derived *>(this);
    }
    catch (exceptions::add_arguments &e)
    {
      std::cout << "Exception caught, returning self in generic merge series." << std::endl;
      return *static_cast<Derived *>(this);
    }
  }

  /// Add generic entity.
  /**
   * Creates series from entity, then adds it.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::add_generic(const T &x)
  {
    // FIXME: replace this ctor (as well as below) with direct ctor from x?
    Derived tmp=Derived(cf_type(x),*static_cast<Derived *>(this));
    return add(tmp);
  }

  /// Subtract generic entity.
  /**
   * Creates series from entity, then subtracts it.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::subtract_generic(const T &x)
  {
    Derived tmp=Derived(cf_type(x),*static_cast<Derived *>(this));
    return subtract(tmp);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Cf2, class LightTermPair>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::term_by_term_multiplication_trig(
    const term_type &t1, const Term<Cf2,trig_type,allocator_type> &t2, LightTermPair &term_pair,
    cf_type &new_c)
  {
    t1.trig().trigmult(t2.trig(),term_pair.template get
      <0>().trig(),
      term_pair.template get<1>().trig());
    term_pair.template get
      <0>().cf()=term_pair.template get
      <1>().cf()=new_c;
    if (t1.trig().flavour() == t2.trig().flavour())
    {
      term_pair.template get
        <0>().trig().flavour()=term_pair.template get
        <1>().trig().flavour()=true;
      if(!t1.trig().flavour())
      {
        term_pair.template get
          <1>().cf().invert_sign();
      }
    }
    else
    {
      term_pair.template get
        <0>().trig().flavour()=term_pair.template get
        <1>().trig().flavour()=false;
      if(t1.trig().flavour())
      {
        term_pair.template get
          <0>().cf().invert_sign();
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
    if (!math::is_zero_vec(lin_args()) or !math::is_zero_vec(ps2.lin_args()))
    {
      std::cout << "Non-zero linargs in multiplication." << std::endl;
      std::exit(1);
      return false;
    }
    try
      { merge_args(ps2);}
    catch (exceptions::add_arguments &e)
    {
      std::cout << "Unable to merge arguments in series_multiplication_preliminaries." << std::endl;
      return false;
    }
    retval.merge_args(*static_cast<Derived *>(this));
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
      static_cast<Derived *>(this)->cf_multiplication(ps2.g_s_index().begin()->cf());
      return true;
    }
    else if (is_cf())
    {
      cf_type tmp(g_s_index().begin()->cf());
      assign_series(ps2);
      static_cast<Derived *>(this)->cf_multiplication(tmp);
      std::cout << "Cf2\n";
      return true;
    }
    return false;
  }

  /// Multiplication by a generic series.
  /**
   * Requires some methods to be implemented in the derived classes.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::mult_by_series(const Derived2 &ps2)
  {
    Derived *derived_cast=static_cast<Derived *>(this);
    Derived retval;
    if (series_multiplication_preliminaries(ps2,retval))
    {
      if (series_multiplication_optimize_for_cf(ps2))
      {
        return *static_cast<Derived *>(this);
      }
      derived_cast->multiply_terms(ps2,retval);
    }
    swap(retval);
    return *static_cast<Derived *>(this);
  }

  /// Multiplication by a generic entity.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::mult_by_generic(const T &c)
  {
    // TODO: maybe this can be share with generic_division with the help of functors?
    if (is_empty())
    {
      return *static_cast<Derived *>(this);
    }
    if (!math::is_zero_vec(lin_args()))
    {
      std::cout << "Non-zero linargs in generic series multiplication." << std::endl;
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
      tmp_term.cf().mult_by(c);
      it_hint=tmp_ps.insert(tmp_term,it_hint);
    }
    swap(tmp_ps);
    return *static_cast<Derived *>(this);
  }

  /// Multiplication by an integer.
  /**
   * This is a bit more complicated because we have to take care of lin_args.
   */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::mult_by(int n)
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
    mult_by_generic(n);
    // Multiply the old linargs and restore them
    for (j=0;j<w;++j)
    {
      lin_args()[j]=old_lin_args[j]*n;
    }
    return *static_cast<Derived *>(this);
  }

  /// Mult by self.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::mult_by(const Derived &ps2)
  {
    return mult_by_series(ps2);
  }

  /// Generic division.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::divide_by_generic(const T &x)
  {
    // NOTICE: can this check be improved?
    if (std::abs(x) == 0)
    {
      std::cout << "ERROR: division by zero in divide_by_generic, returning self." << std::endl;
      std::abort();
      return *static_cast<Derived *>(this);
    }
    if (is_empty())
    {
      return *static_cast<Derived *>(this);
    }
    if (!math::is_zero_vec(lin_args()))
    {
      // NOTICE: maybe here we could deal with exact int/int divisions. Not really important
      // ATM though.
      std::cout << "Non-zero linargs in divide_by_generic!" << std::endl;
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
      tmp_term.cf().divide_by(x);
      it_hint=tmp_ps.template insert<false,true>(tmp_term,it_hint);
    }
    swap(tmp_ps);
    return *static_cast<Derived *>(this);
  }

  /// Generic series addition.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::add_series(const Derived2 &ps2)
  {
    return merge_with_series<Derived2,true>(ps2);
  }

  /// Generic series subtraction.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class Derived2>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::subtract_series(const Derived2 &ps2)
  {
    return merge_with_series<Derived2,false>(ps2);
  }

  /// Add self.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::add(const Derived &ps2)
  {
    return add_series(ps2);
  }

  /// Add int.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::add(int n)
  {
    return add_generic(n);
  }

  /// Add double.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::add(const double &x)
  {
    return add_generic(x);
  }

  /// Subtract int.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::subtract(int n)
  {
    return subtract_generic(n);
  }

  /// Subtract double.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::subtract(const double &x)
  {
    return subtract_generic(x);
  }

  /// Subtract self.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::subtract(const Derived &ps2)
  {
    return subtract_series(ps2);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::mult_by(const double &x)
  {
    return mult_by_generic(x);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::divide_by(int n)
  {
    return divide_by_generic(n);
  }

  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived &base_pseries<__PIRANHA_BASE_PS_TP>::divide_by(const double &x)
  {
    return divide_by_generic(x);
  }
}
#endif
