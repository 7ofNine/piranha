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

#ifndef PIRANHA_BASE_POLYNOMIAL_H
#define PIRANHA_BASE_POLYNOMIAL_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <gmp.h>
#include <gmpxx.h>

#include "math.h"
#include "monomial_gmp_array.h"
#include "symbol_limiter.h"

namespace piranha
{
/// Base polynomial class.
/**
 * Implemented as a sorted and hashed multiindex container for monomials M.
 */
  template <class T, class Derived>
    class base_polynomial
  {
    public:
      typedef monomial_gmp_array<T> m_type;
      typedef typename m_type::expo_type expo_type;
      typedef typename m_type::degree_type degree_type;
    private:
/// Tag for the degree index.
      struct degree
        {}
      ;
/// Tag for the hashed index.
      struct hash
        {}
      ;
      typedef boost::multi_index::indexed_by <
        boost::multi_index::ordered_non_unique <
        boost::multi_index::tag<degree>,
        boost::multi_index::composite_key <
        m_type,
        boost::multi_index::const_mem_fun < m_type, expo_type,
        &m_type::g_min_expo > ,
        boost::multi_index::const_mem_fun < m_type, degree_type,
        &m_type::g_degree > >
        >,
        boost::multi_index::hashed_unique <
        boost::multi_index::tag<hash>,
        boost::multi_index::identity<m_type> >
        > index_type;
    public:
      typedef boost::multi_index_container <m_type, index_type> set_type;
      typedef typename set_type::template index<degree>
        ::type degree_index;
      typedef typename degree_index::iterator it_d_index;
      typedef typename set_type::template index<hash>
        ::type hashed_index;
      typedef typename hashed_index::iterator it_h_index;
      friend class std::complex<Derived>;
// Start INTERFACE definition.
//-------------------------------------------------------
/// Alias for the evaluation type.
      typedef typename m_type::eval_type eval_type;
// Ctors.
/// Default constructor.
      explicit base_polynomial():set_()
        {}
/// Constructor from int.
      explicit base_polynomial(int n)
      {
        insert(m_type(n));
      }
/// Constructor from double.
      explicit base_polynomial(const double &x)
      {
        insert(m_type(x));
      }
/// Constructor from psymbol.
      explicit base_polynomial(const psymbol &)
      {
        m_type tmp_m((size_t)1);
        tmp_m.g_rational_cf()=1;
// TODO: we can ditch the extra ctor here if we allow for assignment of complex double_cf from real.
        tmp_m.g_numerical_cf()=typename m_type::numerical_type(1.);
        tmp_m.container()[0]=1;
        insert(tmp_m);
      }
/// Copy constructor.
      base_polynomial(const base_polynomial &p):set_(p.set_)
        {}
      explicit base_polynomial(const std::string &);
/// Destructor.
      ~base_polynomial()
        {}
// I/O.
      void print_plain(std::ostream &, const vector_psym_p &) const;
      void print_latex(std::ostream &, const vector_psym_p &) const;
// Manipulation.
      void increase_size(const size_t &);
      void append_args(const size_t &);
      void prepend_args(const size_t &);
      void swap(Derived &);
      void clear()
      {
        set_.clear();
      }
// Evaluation
      eval_type t_eval(const double &, const vector_psym_p &) const;
// Probing
      size_t length() const
      {
        return set_.size();
      }
      bool checkup(const size_t &) const;
      bool is_zero(const vector_psym_p &) const;
      double norm(const vector_psym_p &) const;
      degree_type g_degree() const
      {
        if (set_.empty())
        {
          return 0;
        }
        else
        {
          return d_index().begin()->g_degree();
        }
      }
      expo_type g_min_expo() const
      {
        if (set_.empty())
        {
          return 0;
        }
        else
        {
          return d_index().begin()->g_min_expo();
        }
      }
/// Check whether base_polynomial is larger than size w.
      bool larger(const size_t &w) const
      {
        return (width()>w);
      }
/// Check whether base_polynomial is smaller than size w.
      bool smaller(const size_t &w) const
      {
        return (width()<w);
      }
/// Check whether base_polynomial is compatible with size w.
      bool compatible(const size_t &w) const
      {
        return (width()==w);
      }
/// Get actual size of base_polynomial.
      size_t actual_width() const
      {
        return width();
      }
// Maths.
/// Partial derivative with respect to nth argument.
/**
 * Result is stored in input parameter "retval". No assumptions on retval are made.
 * @param[out] retval, Derived class in which the result will be stored.
 */
      void partial(const size_t &n, Derived &retval) const
      {
        retval.clear();
        m_type tmp_m;
        const size_t w=width();
        tmp_m.reserve(w);
        const it_d_index it_f=d_index().end();
        for (it_d_index it=d_index().begin();it!=it_f;++it)
        {
          it->partial(n,tmp_m);
          retval.insert(tmp_m);
        }
      }
// End INTERFACE definition.
//-------------------------------------------------------
      const degree_index &d_index() const
      {
        return set_.template get
          <degree>();
      }
      const hashed_index &h_index() const
      {
        return set_.template get
          <hash>();
      }
      bool empty() const
      {
        return set_.empty();
      }
      template <class Derived2>
        void mult_by_self(const Derived2 &, const vec_expo_index_limit &);
    protected:
      void insert(const m_type &, bool sign=true);
      degree_index &d_index()
      {
        return set_.template get
          <degree>();
      }
      hashed_index &h_index()
      {
        return set_.template get
          <hash>();
      }
      size_t width() const
      {
        if (empty())
        {
          return 0;
        }
        return d_index().begin()->width();
      }
/// Check if the polynomial can be raised to real power.
      bool pow_preliminary_checks() const
      {
// Requisite: first monomial of the polynomial must be purely numerical and positive.
        if (length() > 0)
        {
          if (d_index().begin()->is_symbolic() || d_index().begin()->g_numerical_cf().value() < 0)
          {
            std::cout << "Error: in raising polynomial to power, first monomial must be numerical and positive."
              << std::endl;
            std::abort();
            return false;
          }
        }
        return true;
      }
      bool pow_handle_special_cases(const double &x)
      {
// Special handling in case of empty polynomial.
        if (empty())
        {
          if (std::abs(x)<settings_manager::numerical_zero())
          {
            std::cout << "ERROR: want to calculate 0^0 in polynomial." << std::endl;
            std::abort();
            return true;
          }
          if (x<0)
          {
            std::cout << "ERROR: negative power of zero in polynomial." << std::endl;
            std::abort();
            return true;
          }
// 0^x, just return.
          return true;
        }
        if (std::abs(1-x)<settings_manager::numerical_zero())
        {
// this^1.
          return true;
        }
        if (std::abs(x)<settings_manager::numerical_zero())
        {
// this^0.
          clear();
          insert(m_type((int)1));
          return true;
        }
// Optimization: if the first monomial, which is numerical, is also the only one, simply compute power.
        if (length()==1)
        {
          Derived retval;
          retval.insert(m_type(std::pow(d_index().begin()->g_numerical_cf()*d_index().begin()->g_rational_cf().get_d(),x)));
          swap(retval);
          return true;
        }
        return false;
      }
/*      void pow_suitable(boost::tuple<bool,int> &s, const size_t &n) const
      {
        p_assert(n<width());
// If assert below fails, the retval tuple may end up uninitialized.
        p_assert(length()>=1);
// The minimum index until now is the first monomials's.
        s.get<1>()=d_index().begin()->container()[n];
        int candidate;
        const it_d_index it_f=d_index().end();
        for (it_d_index it=d_index().begin();it!=it_f;++it)
        {
          candidate=it->container()[n];
          if (candidate<=0)
          {
            s.get<0>()=false;
            std::cout << "Non suitability for real powers detected." << std::endl;
            return;
          }
          else
          {
            if (candidate < s.get<1>())
            {
              s.get<1>()=candidate;
            }
          }
        }
        s.get<0>()=true;
        std::cout << "Final minimum index in real power of polynomial is: " << s.get<1>() << std::endl;
      }
*/
/// Maximum natural power allowed.
/**
 * Maximum natural power that will produce a non-null single-variable polynomial,
 * given the desired exponent limit for the variable and the variable's minimum exponent
 * in the polynomial.
 */
/*      unsigned int max_natural_power(int limit, int min_expo) const
      {
        p_assert(limit >= 0 && min_expo > 0 );
        return (unsigned int)std::floor(((float)limit)/min_expo);
      }
*/
      void base_merge(const base_polynomial &, bool);
      template <class N>
        void ll_generic_integer_division(const N &);
      void basic_assignment(const base_polynomial &p)
      {
        if (&p!=this)
        {
          set_=p.set_;
        }
      }
// TODO check if these are really used.
      void mult_by_int(int);
      void mult_by_double(const double &);
      template <class Derived2>
        void mult_by_self(const Derived2 &);
      //void basic_pow(const double &, const vec_expo_index_limit limits &);
    protected:
      static const std::string  separator_;
      set_type                  set_;
  }
  ;

  template <class T, class Derived>
    const std::string base_polynomial<T,Derived>::separator_="&";

/// Constructor from string.
  template <class T, class Derived>
    inline base_polynomial<T,Derived>::base_polynomial(const std::string &str)
  {
    deque_string split_v;
    boost::split(split_v,str,boost::is_any_of(separator_));
    const size_t w=split_v.size();
    for (size_t i=0;i<w;++i)
    {
      boost::trim_left_if(split_v[i],boost::is_any_of("{"));
      boost::trim_right_if(split_v[i],boost::is_any_of("}"));
      m_type tmp_m(split_v[i]);
      if (tmp_m.larger(width()) && !empty())
      {
        increase_size(tmp_m.width());
      }
      insert(tmp_m);
    }
  }

/// Print in plain format.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::print_plain(std::ostream &out_stream, const vector_psym_p &cv) const
  {
// TODO: cache end() here.
    stream_manager::setup_print(out_stream);
    out_stream << '{';
    for (it_d_index it=d_index().begin();it!=d_index().end();++it)
    {
      it->print_plain(out_stream,cv);
      ++it;
      if (it!=d_index().end())
      {
        out_stream << '&';
      }
      --it;
    }
    out_stream << '}';
  }

/// Print in latex format.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::print_latex(std::ostream &out_stream, const vector_psym_p &cv) const
  {
// TODO: cache end() here.
    stream_manager::setup_print(out_stream);
    for (it_d_index it=d_index().begin();it!=d_index().end();++it)
    {
      if (it->g_numerical_cf().is_positive() && it!=d_index().begin())
      {
        out_stream << std::string("$+$");
      }
      it->print_latex(out_stream,cv);
    }
  }

/// Increase size.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::increase_size(const size_t &w)
  {
    p_assert(w>=width());
    Derived retval;
    const it_d_index it_f=d_index().end();
    for (it_d_index it=d_index().begin();it!=it_f;++it)
    {
      m_type tmp_m(*it);
      tmp_m.increase_size(w);
      retval.insert(tmp_m);
    }
    swap(retval);
  }

/// Append arguments.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::append_args(const size_t &n)
  {
    increase_size(width()+n);
  }

/// Prepend arguments.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::prepend_args(const size_t &n)
  {
    Derived retval;
    const it_d_index it_f=d_index().end();
    for (it_d_index it=d_index().begin();it!=it_f;++it)
    {
      m_type tmp_m(*it);
      tmp_m.prepend_args(n);
      retval.insert(tmp_m);
    }
    swap(retval);
  }

/// Swap contents with another base_polynomial.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::swap(Derived &p)
  {
    set_.swap(p.set_);
  }

/// Insert a new element.
  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::insert(const m_type &m, bool sign)
  {
    if (m.is_zero())
    {
      return;
    }
    const size_t w=width();
    m_type *new_m=0;
    if (m.smaller(w))
    {
// FIXME: used boost shared pointer.
      new_m = new m_type(m);
      new_m->increase_size(w);
    }
    p_assert(!(m.larger(w) && !empty()))
/*
    if (!set_.empty() && m.larger(w))
    {
      increase_size(m.width());
    }
*/
    const m_type *insert_m;
    if (new_m==0)
    {
      insert_m=&m;
    }
    else
    {
      insert_m=new_m;
    }
    it_h_index it=h_index().find(*insert_m);
    if (it==h_index().end())
    {
// The term is NOT a duplicate, insert in the set.
      std::pair<it_h_index,bool> result=h_index().insert(*insert_m);
      p_assert(result.second);
// If requested, subtract.
      if (!sign)
      {
        p_assert(h_index().modify(result.first,
          typename m_type::update_numerical_cf(-result.first->g_numerical_cf())));
      }
    }
    else
    {
// The term is in the set, hence an existing term will be modified.
      typename m_type::numerical_type numerical_cf;
      typename m_type::rational_type rational_cf;
// Add or subtract according to request.
      m_type::cfs_algebraic_sum(
        it->g_numerical_cf(),insert_m->g_numerical_cf(),it->g_rational_cf(),insert_m->g_rational_cf(),numerical_cf,
        rational_cf,sign);
// Check if the resulting coefficient can be ignored (ie it is small).
      if (numerical_cf.abs()<settings_manager::numerical_zero() || rational_cf==0)
      {
        h_index().erase(it);
      }
      else
      {
        p_assert(h_index().modify(it,typename m_type::update_cfs(numerical_cf,rational_cf)));
      }
    }
    delete new_m;
  }

  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::base_merge(const base_polynomial &p, bool op)
  {
    if (&p==this)
    {
      std::cout << "Need to duplicate!" << std::endl;
      base_merge(base_polynomial(p),op);
    }
    else
    {
      it_h_index it_f=p.h_index().end();
      for (it_h_index it=p.h_index().begin();it!=it_f;++it)
      {
        insert(*it,op);
      }
    }
  }

/// Multiply by self, simple implementation.
/**
 * Multiply all monomials without truncations.
 */
  template <class T, class Derived>
    template <class Derived2>
    inline void base_polynomial<T,Derived>::mult_by_self(const Derived2 &p)
  {
    if ((void *)&p==(void *)this)
    {
      mult_by_self(Derived2(p));
      return;
    }
    Derived retval;
    m_type temp_m;
    const it_h_index it_f1=h_index().end();
    const typename Derived2::it_h_index it_f2=p.h_index().end();
    typename Derived2::it_h_index it2;
    for (it_h_index it1=h_index().begin();it1!=it_f1;++it1)
    {
      for (it2=p.h_index().begin();it2!=it_f2;++it2)
      {
        it1->mult_by(*it2,temp_m);
        retval.insert(temp_m);
      }
    }
    swap(retval);
  }

/// Multiply by self limiting the exponents of symbols.
/**
 * Exponents limits are usually fetched from piranha::symbol_limiter.
 * @see piranha::symbol_limiter.
 */
  template <class T, class Derived>
    template <class Derived2>
    inline void base_polynomial<T,Derived>::mult_by_self(const Derived2 &p,
    const vec_expo_index_limit &v)
  {
// This function is to be used only from Poisson series, in which merging of arguments should
// ensure the validity of the following assert.
    p_assert(width()>=p.width());
    if ((void *)&p==(void *)this)
    {
      mult_by_self(Derived2(p),v);
      return;
    }
    Derived retval;
    m_type temp_m;
    const it_h_index it_f1=h_index().end();
    const typename Derived2::it_h_index it_f2=p.h_index().end();
    typename Derived2::it_h_index it2;
    const size_t w=v.size();
    p_assert(w<=width());
    size_t j;
    expo_type ex1, ex2;
    bool proceed;
    for (it_h_index it1=h_index().begin();it1!=it_f1;++it1)
    {
      for (it2=p.h_index().begin();it2!=it_f2;++it2)
      {
        proceed=true;
        for (j=0;j<w;++j)
        {
// Find the exponents of the limited arguments.
          const size_t index=v[j].get<0>();
// Make sure we are not going outside boundaries.
          p_assert(index < width());
          ex1=it1->g_container()[index];
// it2 can be smaller, since we support multiplication by a narrower polynomial.
          if (it2->smaller(index))
          {
            ex2=0;
          }
          else
          {
            ex2=it2->g_container()[index];
          }
          if ((ex1+ex2)>v[j].get<1>())
          {
            proceed=false;
            break;
          }
        }
        if (proceed)
        {
          it1->mult_by(*it2,temp_m);
          retval.insert(temp_m);
        }
      }
    }
    swap(retval);
  }

  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::mult_by_int(int n)
  {
    if (n==0)
    {
      set_.clear();
      return;
    }
    const it_h_index it_f=h_index().end();
// TODO: Move "if" outside "for" and make two cycles to optimize?
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
// Distinguish between positive and negative, we want only the numerical coefficient to be signed.
      if (n>0)
      {
        p_assert(h_index().modify(it,typename m_type::update_rational_cf(n*it->g_rational_cf())));
      }
      else
      {
        p_assert(h_index().modify(it,typename m_type::update_cfs(-it->g_numerical_cf(),
          (-n)*it->g_rational_cf())));
      }
    }
  }

  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::mult_by_double(const double &x)
  {
    if (std::abs(x)<settings_manager::numerical_zero())
    {
      set_.clear();
      return;
    }
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      p_assert(h_index().modify(it,typename m_type::update_numerical_cf(it->g_numerical_cf()*x)));
    }
  }

// Low-level generic division by integer type.
  template <class T, class Derived>
    template <class N>
    inline void base_polynomial<T,Derived>::ll_generic_integer_division(const N &n)
  {
    if (n==0)
    {
      std::cout << "FATAL: division by zero." << std::endl;
      std::exit(1);
    }
    const it_h_index it_f=h_index().end();
// TODO: Move "if" outside "for" and make two cycles to optimize?
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
// Distinguish between positive and negative, we want only the numerical coefficient to be signed.
      if (n>0)
      {
        p_assert(h_index().modify(it,typename m_type::update_rational_cf(it->g_rational_cf()/n)));
      }
      else
      {
        p_assert(h_index().modify(it,typename m_type::update_cfs(-it->g_numerical_cf(),
          it->g_rational_cf()/(-n))));
      }
    }
  }

/// Time evaluation.
  template <class T, class Derived>
    inline typename base_polynomial<T,Derived>::eval_type base_polynomial<T,Derived>::t_eval(const double &t, const vector_psym_p &v) const
  {
    p_assert(width()==v.size());
    eval_type retval=0;
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      retval+=it->t_eval(t,v);
    }
    return retval;
  }

/// Diagnostic checkup.
  template <class T, class Derived>
    inline bool base_polynomial<T,Derived>::checkup(const size_t &s) const
  {
    const it_h_index it_f=h_index().end();
    size_t w=width();
// Let's check that base_polynomial's size is equal to s.
    if (s!=w)
    {
      std::cout << "Size mismatch in base_polynomial." << std::endl;
    }
// Let's check that all monomials have same size.
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      if (!it->checkup(w))
      {
        std::cout << "Size mismatch in monomial." << std::endl;
        return false;
      }
    }
    return true;
  }

/// Check whether a base_polynomial is zero.
// FIXME: do we really need vector of symbols as input here?
  template <class T, class Derived>
    inline bool base_polynomial<T,Derived>::is_zero(const vector_psym_p &v) const
  {
    if (empty())
    {
      return true;
    }
// We have to check that monomials do not all evaluate to zero... In case there are
// symbols with zero constant part.
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      if (std::abs(it->t_eval(0,v))>settings_manager::numerical_zero())
      {
        return false;
      }
    }
    return true;
  }

/// Get base_polynomial's norm.
  template <class T, class Derived>
    inline double base_polynomial<T,Derived>::norm(const vector_psym_p &v) const
  {
    eval_type retval=0.;
    const it_h_index it_f=h_index().end();
    for (it_h_index it=h_index().begin();it!=it_f;++it)
    {
      retval+=it->t_eval(0,v);
    }
    return std::abs(retval);
  }

/// Basic real power.
/*  template <class T, class Derived>
    inline void base_polynomial<T,Derived>::basic_pow(const double &x, const vec_expo_index_limit limits &v)
  {
    if (!pow_preliminary_checks())
    {
      return;
    }
    if (pow_handle_special_cases(x))
    {
      return;
    }
// This is the remainder part of this, i.e. this without the leading term.
    Derived x(*static_cast<Derived *>(this));
    x.d_index().erase(x.d_index().begin());
// Now we have to check if with the provided vector for symbol limits we effectively reach a point after which
// the binomial expansion will lead to null polynomials being created. We must stop at that point.
// In binomial expansions, "x" is raised to a power equal to the number of iterations.
    int iterations=-1;
    const size_t v_size=v.size();
    boost::tuple<bool,int> is_suitable(true,0);
    for (size_t i=0;j<v_size;++j)
    {
// Make sure we are not poking outside the boundaries of the monomial's container for exponents.
      p_assert(v[i].get<0>()<x.width());
      x.pow_suitable(is_suitable);
      if (!is_suitable.get<0>())
      {
        std::cout << "Polynomial cannot be raised to real power given current exponents limits, returning self." << std::endl;
        std::abort();
        return;
      }
      if (v[i].get<1>() < 0)
      {
        std::cout << "Cannot accept a negative symbol limit when raising to power." << std::endl;
        std::abort();
        return;
      }
max_natural_power(v[i].get<1>(),)
      if (is_suitable.get<1>() > iterations)
      {
        iterations=is_suitable.get<1>();
      }
    }
    if (iterations < 0)
    {
      std::cout << "Negative iterations for power, returning self;
      std::abort();
      return;
    }


max_natural_power(int limit, int min_expo) const



// NOTICE: Hard coded binomial expansion error to 1/10 of desired precision.
    const double error=.1*std::pow(derived_cast->g_norm(),power)*
      settings_manager::prec();
    const unsigned int limit_index=pow_limit(error,power);
    Derived retval, x(*derived_cast), tmp(1.);
    x.term_erase(x.s_index().begin());
    for (unsigned int i=0;i<=limit_index;++i)
    {
      {
        Derived tmp2(tmp);
        tmp2*=math::choose(power,i);
        tmp2*=Derived(a.pow(power-i),tmp2);
        retval+=tmp2;
      }
      tmp*=x;
    }
    return retval;














    if (set_.size()>1)
    {
      std::cout << "ERROR: won't invert a non singular polynomial." << std::endl;
      std::abort();
    }
    Derived retval;
// Real case.
    retval.insert(d_index().begin()->pow(x));
    swap(retval);
  }*/
}
#endif
