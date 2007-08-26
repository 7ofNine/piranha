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

#ifndef PIRANHA_MONOMIAL_GMP_ARRAY_H
#define PIRANHA_MONOMIAL_GMP_ARRAY_H

#include <boost/functional/hash/hash.hpp>         // hash_combine.
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <sstream>

#include "common_typedefs.h"
#include "p_assert.h"
#include "psymbol.h"
#include "utils.h"
#include "math.h"

namespace piranha
{
// FIXME: fix the order of data members (small ones first!)
/// Base monomial class, array implementation.
/**
 * Monomial is made of three parts: a numerical coefficient, a rational coefficient implemented
 * with the GMP mpq_class and an array for the exponents of the variables.
 */
  template <class T>
    class monomial_gmp_array
  {
    public:
// Start INTERFACE definition.
//-------------------------------------------------------
/// Alias for the numerical coefficient.
      typedef T numerical_type;
/// Alias for the rational coefficient.
      typedef mpq_class rational_type;
/// Alias for the container type.
      typedef std::valarray<expo_type> container_type;
/// Evaluation result.
      typedef typename numerical_type::eval_type eval_type;
/// Functor to update monomial's coefficients.
      struct update_cfs
      {
        update_cfs(const numerical_type &numerical, const rational_type &rational):
        numerical_(numerical),rational_(rational)
          {}
        ~update_cfs()
          {}
        void operator()(monomial_gmp_array &m) const
        {
          m.private_numerical_cf_=numerical_;
          m.private_rational_cf_=rational_;
        }
        const numerical_type  &numerical_;
        const rational_type   &rational_;
      }
      ;
/// Functor to update monomial's numerical coefficient.
      struct update_numerical_cf
      {
        update_numerical_cf(const numerical_type &numerical):
        numerical_(numerical)
          {}
        ~update_numerical_cf()
          {}
        void operator()(monomial_gmp_array &m) const
        {
          m.private_numerical_cf_=numerical_;
        }
        const numerical_type   &numerical_;
      }
      ;
/// Functor to update monomial's rational coefficient.
      struct update_rational_cf
      {
        update_rational_cf(const rational_type &rational):
        rational_(rational)
          {}
        ~update_rational_cf()
          {}
        void operator()(monomial_gmp_array &m) const
        {
          m.private_rational_cf_=rational_;
        }
        const rational_type   &rational_;
      }
      ;
// Ctors.
/// Default constructor.
      explicit monomial_gmp_array():private_min_expo_(0),private_degree_(0),private_container_(0),
        private_numerical_cf_(0.),private_rational_cf_(0)
        {}
/// Constructor from size.
      explicit monomial_gmp_array(const size_t &n):private_min_expo_(0),private_degree_(0),private_container_(n),
        private_numerical_cf_(0.),private_rational_cf_(0)
        {}
/// Constructor from integer.
      explicit monomial_gmp_array(int n):private_min_expo_(0),private_degree_(0),private_container_(0),
        private_numerical_cf_(math::sgn(n)),private_rational_cf_((unsigned int)std::abs(n))
      {
        if (n==0)
        {
          private_numerical_cf_=numerical_type(0.);
          private_rational_cf_=1;
        }
      }
/// Constructor from numerical type.
      explicit monomial_gmp_array(const numerical_type &x):private_min_expo_(0),private_degree_(0),private_container_(0),
        private_numerical_cf_(x),private_rational_cf_(1)
        {}
/// Constructor from double.
      explicit monomial_gmp_array(const double &x):private_min_expo_(0),private_degree_(0),private_container_(0),
        private_numerical_cf_(x),private_rational_cf_(1)
        {}
      explicit monomial_gmp_array(const std::string &);
/*/// Constructor from numerical, rational and exponent container.
      explicit monomial_gmp_array(const numerical_type &num, const rational_type &rat, const container_type &cont)
        :private_container_(cont),private_numerical_cf_(num),private_rational_cf_(rat)
        {}*/
/// Copy constructor.
      monomial_gmp_array(const monomial_gmp_array &m):private_min_expo_(m.private_min_expo_),private_degree_(m.private_degree_),
        private_container_(m.private_container_),private_numerical_cf_(m.private_numerical_cf_),
        private_rational_cf_(m.private_rational_cf_)
        {}
/// Copy constructor from monomial with different numerical type.
      template <class U>
        monomial_gmp_array(const monomial_gmp_array<U> &m):private_min_expo_(m.g_min_expo()),private_degree_(m.get_degree()),
        private_container_(m.g_container()),private_numerical_cf_(m.g_numerical_cf()),private_rational_cf_(m.g_rational_cf())
        {}
/// Constructor from psymbol.
      monomial_gmp_array(const psymbol &):private_min_expo_(1),private_degree_(1),private_container_(size_t(1)),
        private_numerical_cf_(1.),private_rational_cf_(1)
      {
        private_container_[0]=1;
      }
/// Destructor.
      ~monomial_gmp_array()
        {}
// Getters & Setters.
/// Get minimum exponent.
      expo_type g_min_expo() const
      {
        return private_min_expo_;
      }
/// Get degree.
      degree_type g_degree() const
      {
        return private_degree_;
      }
/// Get size of monomial (i.e., the number of variables).
      size_t width() const
      {
        return g_container().size();
      }
/// Get const reference to the numerical coefficient.
      const numerical_type &g_numerical_cf() const
      {
        return private_numerical_cf_;
      }
/// Get reference to the numerical coefficient.
      numerical_type &s_numerical_cf()
      {
        return private_numerical_cf_;
      }
/// Get const reference to the rational coefficient.
      const rational_type &g_rational_cf() const
      {
        return private_rational_cf_;
      }
/// Get reference to the rational coefficient.
      rational_type &s_rational_cf()
      {
        return private_rational_cf_;
      }
/// Get const reference to the container.
      const container_type &g_container() const
      {
        return private_container_;
      }
// I/O.
      void print_plain(std::ostream &, const vector_psym_p &) const;
      void print_latex(std::ostream &, const vector_psym_p &) const;
// Manip.
      void increase_size(const size_t &);
      void prepend_args(const size_t &);
// Evaluation.
      eval_type t_eval(const double &, const vector_psym_p &) const;
// Probing.
/// Explicitly calculate minimum exponent.
      expo_type calc_min_expo() const
      {
        const size_t w=width();
        if (w==0)
        {
          return 0;
        }
        expo_type retval=g_container()[0];
        for (size_t i=1;i<w;++i)
        {
          if (g_container()[i] < retval)
          {
            retval = g_container()[i];
          }
        }
        return retval;
      }
/// Explictly calculate degree.
      degree_type calc_degree() const
      {
        const size_t w=width();
        degree_type retval=0;
        for (size_t i=0;i<w;++i)
        {
          retval+=g_container()[i];
        }
        return retval;
      }
      size_t hasher() const;
      bool smaller(const size_t &) const;
      bool larger(const size_t &) const;
      bool is_zero() const;
      bool checkup(const size_t &) const;
      bool operator==(const monomial_gmp_array &) const;
// Maths.
      template <class U>
        void mult_by(const monomial_gmp_array<U> &, monomial_gmp_array &) const;
      static void cfs_algebraic_sum(const numerical_type &, const numerical_type &,
        const rational_type &, const rational_type &, numerical_type &, rational_type &, bool);
// TODO: these two should provide a retval as parameter.
      monomial_gmp_array pow(int) const;
      monomial_gmp_array pow(const double &) const;
/// Partial derivative.
/**
 * Calculated with respect to argument at position n. Retval is assumed to be able to accept the result,
 * i.e. its exponent container must be the right size.
 * @param[out] retval: piranha::monomial_gmp_array where result is stored.
 */
      void partial(const size_t &n, monomial_gmp_array &retval) const
      {
        const size_t w=width();
        p_assert(n<w);
        retval.increase_size(w);
        p_assert(w==retval.width());
        const expo_type exponent=g_container()[n];
        retval.private_min_expo_=private_min_expo_;
        retval.private_degree_=private_degree_;
        if (exponent==0)
        {
// If the exponent is zero, partial derivative will be zero too.
          retval.s_numerical_cf()=numerical_type(0);
        }
        else
        {
// The rational cf of retval gets multiplied by the exponent.
          retval.s_rational_cf()=g_rational_cf();
// Take abs because only numerical cf can be negative.
          retval.s_rational_cf()*=(expo_type)std::abs(exponent);
          if (exponent>0)
          {
            retval.s_numerical_cf()=g_numerical_cf();
          }
          else
          {
            retval.s_numerical_cf()=-g_numerical_cf();
          }
// Assign same exponents as this...
          retval.private_container_=private_container_;
// ... just modify the interested one.
          --retval.private_container_[n];
// Take care of degree and minimum exponent.
          --retval.private_degree_;
          if (retval.private_container_[n] < retval.private_min_expo_)
          {
            retval.private_min_expo_=retval.private_container_[n];
          }
// TODO: place assert to make sure we don't go out expo_type range?
// This should be generalized, if we decide this way.
        }
      }
// Operators.
      template <class U>
        monomial_gmp_array &operator=(const monomial_gmp_array<U> &);
      bool operator<(const monomial_gmp_array &) const;
// End INTERFACE definition.
//-------------------------------------------------------
    protected:
      template <class U>
        void basic_assignment(const monomial_gmp_array<U> &);
      bool is_symbolic() const;
    private:
      expo_type                   private_min_expo_;
      degree_type                 private_degree_;
      container_type              private_container_;
      numerical_type              private_numerical_cf_;
      rational_type               private_rational_cf_;
// Separator used in I/O.
      static const std::string    separator1_;
      static const std::string    separator2_;
  }
  ;

  template<class T>
    const std::string monomial_gmp_array<T>::separator1_=":";
  template<class T>
    const std::string monomial_gmp_array<T>::separator2_=" ";

/// Constructor from string.
  template <class T>
    inline monomial_gmp_array<T>::monomial_gmp_array(const std::string &s):private_min_expo_(0),private_degree_(0),
    private_container_(0),private_numerical_cf_(0.),private_rational_cf_(0)
  {
    deque_string split_v;
    boost::split(split_v,s,boost::is_any_of(separator1_));
    if (split_v.size()==1)
    {
// In this case we try to cast to double.
      std::cout << "Trying to build monomial from numerical value." << std::endl;
      private_rational_cf_=1;
      private_numerical_cf_=utils::lexical_converter<numerical_type>(split_v[0]);
      return;
    }
    if (split_v.size()!=3)
    {
      std::cout << "Invalid number of elements in monomial_gmp_array string. Returning 0." << std::endl;
      return;
    }
    boost::trim(split_v[0]);
    private_numerical_cf_=utils::lexical_converter<numerical_type>(split_v[0]);
    boost::trim(split_v[1]);
    private_rational_cf_=utils::lexical_converter<rational_type>(split_v[1]);
// TODO: this should be removed if we are going to switch to boost rational class.
    private_rational_cf_.canonicalize();
    if (private_rational_cf_<0)
    {
      private_rational_cf_=abs(private_rational_cf_);
    }
    boost::trim(split_v[2]);
    if (split_v[2].size()<2 || *(split_v[2].begin())!='[' || *(--split_v[2].end())!=']')
    {
      std::cout << "Error converting string to exponents in monomial_gmp_array, returning 0." << std::endl;
      return;
    }
    boost::trim_left_if(split_v[2],boost::is_any_of("["));
    boost::trim_right_if(split_v[2],boost::is_any_of("]"));
    deque_string exponent_v;
    boost::split(exponent_v,split_v[2],boost::is_any_of(separator2_));
    private_container_.resize(exponent_v.size());
    for (size_t i=0;i<exponent_v.size();++i)
    {
      if (exponent_v[i]!="")
      {
        private_container_[i]=utils::lexical_converter<expo_type>(exponent_v[i]);
      }
// TODO: is this necessary or does valarray def-initialize elements to zero?
//       else
//       {
//         private_container_[i]=0;
//       }
// Add newfound exponent to monomial's degree.
      private_degree_+=private_container_[i];
    }
    private_min_expo_=calc_min_expo();
  }

/// Print in plain format.
  template <class T>
    inline void monomial_gmp_array<T>::print_plain(std::ostream &out_stream, const vector_psym_p &) const
  {
    stream_manager::setup_print(out_stream);
    out_stream << g_numerical_cf() << separator1_ << g_rational_cf() << separator1_;
    out_stream << '[';
    for (size_t i=0;i<width();++i)
    {
      out_stream << g_container()[i];
      if (i!=width()-1)
      {
        out_stream << separator2_;
      }
    }
    out_stream << ']';
  }

/// Print in latex format.
  template <class T>
    inline void monomial_gmp_array<T>::print_latex(std::ostream &out_stream, const vector_psym_p &v) const
  {
    p_assert(v.size()==width());
    std::ostringstream tmp;
    stream_manager::setup_print(tmp);
    tmp << std::string("$");
// Print only if it is not 1 in abs and it is not the only element to be printed. If it is -1 print the sign.
    if (g_numerical_cf().is_unity() && (is_symbolic() || g_rational_cf() != 1))
    {
      if (g_numerical_cf().is_negative())
      {
        tmp << std::string("-");
      }
    }
    else
    {
      tmp << g_numerical_cf();
// Print cdot if the following fraction is really an integer - but not 1, because in that
// case it would not be printed.
      if (g_rational_cf().get_den()==1 && g_rational_cf().get_num()!=1)
      {
        tmp << std::string("\\cdot ");
      }
    }
    if (g_rational_cf()!=1)
    {
// If fraction is really an integer print only numerator.
      if (g_rational_cf().get_den()==1)
      {
        tmp << g_rational_cf().get_num();
      }
      else
      {
        tmp << "\\frac{" << g_rational_cf().get_num() << "}{" << g_rational_cf().get_den() << "}";
      }
    }
    for (size_t i=0;i<width();++i)
    {
      if (g_container()[i]==0)
        {}
        else if (g_container()[i]==1)
      {
        tmp << v[i]->name();
      }
      else
      {
// Trailing space needed to avoid problems with symbols whose name starts with a number.
        tmp << v[i]->name() << std::string("^{") << g_container()[i] << "} ";
      }
    }
    tmp << std::string("$");
// Print only if we have something to print.
    if (tmp.str()!=std::string("$$"))
    {
      out_stream << tmp.str();
    }
  }

/// Check whether monomial is symbolic or purely numerical.
  template <class T>
    inline bool monomial_gmp_array<T>::is_symbolic() const
  {
    return (g_min_expo()!=0);
  }

/// Hasher function.
  template <class T>
    inline size_t monomial_gmp_array<T>::hasher() const
  {
    size_t seed=0;
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      boost::hash_combine(seed,g_container()[i]);
    }
    return seed;
  }

/// Overload for the computation of the hash_value.
  template <class T>
    inline size_t hash_value(const monomial_gmp_array<T> &m)
  {
    return m.hasher();
  }

/// Operator less-than.
  template <class T>
    inline bool monomial_gmp_array<T>::operator<(const monomial_gmp_array &m2) const
  {
    p_assert(width()==m2.width());
    const expo_type ex1=g_min_expo(), ex2=m2.g_min_expo();
    if (ex1 < ex2)
    {
      return true;
    }
    if (ex1 > ex2)
    {
      return false;
    }
    const degree_type deg1=g_degree(), deg2=m2.g_degree();
    if (deg1 < deg2)
    {
      return true;
    }
    if (deg1 > deg2)
    {
      return false;
    }
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      if (g_container()[i]<m2.g_container()[i])
      {
        return true;
      }
      if (g_container()[i]>m2.g_container()[i])
      {
        return false;
      }
    }
// The two monomials are equal.
    return false;
  }

/// Basic assignment.
/**
 * Assignation to a smaller monomial will result in assertion failure.
 */
  template <class T>
    template <class U>
    inline void monomial_gmp_array<T>::basic_assignment(const monomial_gmp_array<U> &m2)
  {
    if ((void *)&m2!=(void *)this)
    {
// First we must resize, before assigning: valarray does not do it automatically.
      increase_size(m2.width());
      private_container_=m2.g_container();
      private_numerical_cf_=m2.g_numerical_cf();
      private_rational_cf_=m2.g_rational_cf();
      private_degree_=m2.g_degree();
      private_min_expo_=m2.g_min_expo();
    }
  }

/// Assignment operator.
  template <class T>
    template <class U>
    inline monomial_gmp_array<T> &monomial_gmp_array<T>::
    operator=(const monomial_gmp_array<U> &m2)
  {
    basic_assignment(m2);
    return *this;
  }

/// Check whether monomial is smaller than size w.
  template <class T>
    inline bool monomial_gmp_array<T>::smaller(const size_t &w) const
  {
    return (width()<w);
  }

/// Check whether monomial is larger than size w.
  template <class T>
    inline bool monomial_gmp_array<T>::larger(const size_t &w) const
  {
    return (width()>w);
  }

/// Diagnostic check.
  template <class T>
    inline bool monomial_gmp_array<T>::checkup(const size_t &w) const
  {
    return (width()==w && g_min_expo()==calc_min_expo() && g_degree()==calc_degree());
  }

/// Time evaluation.
  template <class T>
    inline typename monomial_gmp_array<T>::eval_type
    monomial_gmp_array<T>::t_eval(const double &t, const vector_psym_p &v) const
  {
    const size_t w=width();
    p_assert(v.size()==w);
    eval_type retval=g_numerical_cf().value()*g_rational_cf().get_d();
    for (size_t i=0;i<w;++i)
    {
// FIXME: use natural_pow or the like here, to speed up.
      retval*=std::pow(v[i]->t_eval(t),g_container()[i]);
    }
    return retval;
  }

/// Test for equality.
  template <class T>
    inline bool monomial_gmp_array<T>::operator==(const monomial_gmp_array &m2) const
  {
    const size_t w=width();
    p_assert(w==m2.width());
    for (size_t i=0;i<w;++i)
    {
      if (g_container()[i]!=m2.g_container()[i])
      {
        return false;
      }
    }
    p_assert(g_min_expo()==m2.g_min_expo() && g_degree()==m2.g_degree());
    return true;
  }

/// Increrase size so that final size is new_w.
  template <class T>
    inline void monomial_gmp_array<T>::increase_size(const size_t &new_w)
  {
    p_assert(new_w>=width());
    if (new_w>width())
    {
      container_type old=private_container_;
      private_container_.resize(new_w);
      const size_t old_w=old.size();
      for (size_t i=0;i<old_w;++i)
      {
        private_container_[i]=old[i];
      }
// Refresh minimum exponent.
      if (g_min_expo() > 0)
      {
        private_min_expo_=0;
      }
    }
  }

/// Prepend arguments.
  template <class T>
    inline void monomial_gmp_array<T>::prepend_args(const size_t &n)
  {
    if (n>0)
    {
      container_type old=private_container_;
      const size_t old_w=old.size();
      private_container_.resize(old_w+n);
      for (size_t i=0;i<old_w;++i)
      {
        private_container_[i+n]=old[i];
      }
    }
  }

/// Check whether monomial is zero.
  template <class T>
    inline bool monomial_gmp_array<T>::is_zero() const
  {
// We must do this way because if we get_d() immediately on rational_cf_ we may obtain zero for very
// small fractions, but the numerical_cf_ could counter-balance the smallness hence leading to a
// "non-zero" monomial.
// FIXME: here we are copying a gmp rational, it is costly. Check if it is really needed.
    rational_type tmp(g_rational_cf());
    tmp*=std::abs(g_numerical_cf().value());
    p_assert(tmp>=0);
    if (std::abs(tmp.get_d())<settings_manager::numerical_zero())
    {
      return true;
    }
    return false;
  }

/// Multiplication.
/**
 * Multiply by another monomial_gmp_array, and place result in out_m. No assumptions are made on
 * out_m.
 * @param[in] m2 monomial_gmp_array factor.
 * @param[out] out_m monomial_gmp_array that stores the result.
 */
  template <class T>
    template <class U>
    inline void monomial_gmp_array<T>::mult_by(const monomial_gmp_array<U> &m2, monomial_gmp_array &out_m) const
  {
    p_assert(width() >= m2.width());
    out_m.private_numerical_cf_=g_numerical_cf();
    out_m.private_numerical_cf_.mult_by_self(m2.g_numerical_cf());
    out_m.private_rational_cf_=g_rational_cf();
    out_m.private_rational_cf_*=m2.g_rational_cf();
    const size_t w1=width(), w2=m2.width();
    out_m.increase_size(w1);
    size_t i;
    for (i=0;i<w2;++i)
    {
      out_m.private_container_[i]=g_container()[i]+m2.g_container()[i];
    }
    for(;i<w1;++i)
    {
      out_m.private_container_[i]=g_container()[i];
    }
    out_m.private_degree_=g_degree()+m2.g_degree();
    out_m.private_min_expo_=g_min_expo()+m2.g_min_expo();
//std::cout << "OUT_M is: " << out_m << std::endl;
  }

/// Algebraic sum.
  template <class T>
    inline void monomial_gmp_array<T>::cfs_algebraic_sum(const numerical_type &r1, const numerical_type &r2,
    const rational_type &q1, const rational_type &q2, numerical_type &r_ret, rational_type &q_ret, bool sign)
  {
    r_ret=r1;
    q_ret=q1;
// Check whether the numerical parts are "equal".
    if (r1.equal_or_opposite(r2))
    {
// We must check if signs are the same or not.
      bool same_sign=r1.same_sign(r2);
// If (signs are equal and summation) OR (signs are different and subtraction)
      if (same_sign==sign)
      {
        q_ret+=q2;
      }
// If (signs are different and summation) OR (signs are equal and subtraction)
      else
      {
        q_ret-=q2;
// We do not want rationals with signs.
        if (q_ret<0)
        {
          q_ret=abs(q_ret);
        }
      }
    }
    else if (q1==q2)
    {
      if (sign)
      {
        r_ret+=r2;
      }
      else
      {
        r_ret-=r2;
      }
    }
    else
    {
// We do rational division to avoid having a pole later when converting rational to double.
      rational_type tmpq=q2;
      tmpq/=q1;
      if (sign)
      {
        r_ret+=r2*(tmpq.get_d());
      }
      else
      {
        r_ret-=r2*(tmpq.get_d());
      }
    }
  }

/// Integer power.
// TODO: uniform naming of power? i.e., natural power and pow, but here IIRC we need also negative int powers?
// Mmhmm... check out.
// Also: use return values as parameters for inv, powers, etc.
/*  template <class T>
    inline monomial_gmp_array<T> monomial_gmp_array<T>::pow(int n_) const
  {
    p_assert(g_numerical_cf().norm()>settings_manager::numerical_zero());
    p_assert(rational_cf_!=0);
    const size_t w=width();
    monomial_gmp_array retval(w);
    unsigned int n;
    if (n_<0)
    {
      n=(unsigned int)(-n_);
retval.numerical_cf_=math::natural_pow(n,numerical_cf_.inv());
retval.rational_cf_=math::natural_pow(n,
rational_type(rational_cf_.get_den(),rational_cf_.get_num()));
}
else
{
n=(unsigned int)(n_);
retval.numerical_cf_=math::natural_pow(n,numerical_cf_);
retval.rational_cf_=math::natural_pow(n,rational_cf_);
}
// Take care of arguments.
for (size_t i=0;i<w;++i)
{
retval.container_[i]=container_[i]*n_;
}
return retval;
}

/// Real power.
template <class T>
inline monomial_gmp_array<T> monomial_gmp_array<T>::pow(const double &x) const
{
if (is_symbolic())
{
std::cout << "ERROR: real power of symbolic monomial." << std::endl;
std::abort();
}
if (numerical_cf_<0)
{
std::cout << "ERROR: real power of negative numerical coefficient." << std::endl;
std::abort();
}
// NOTICE: here we are counting on the fact that container_ is zero-initialized.
monomial_gmp_array retval(width());
retval.numerical_cf_=(numerical_cf_.pow(x)*=std::pow(rational_cf_.get_d(),x));
retval.rational_cf_=1;
return retval;
}*/
}
#endif
