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
/// Base monomial class, array implementation.
/**
 * Monomial is made of three parts: a numerical coefficient, a rational coefficient implemented
 * with the GMP mpq_class and an array for the exponents of the variables.
 */
  template <class T>
    class monomial_gmp_array
  {
    public:
      typedef short int expo_t;
// Start INTERFACE definition.
//-------------------------------------------------------
// TODO: use "self" type here (and polynomial too)?
/// Alias for the degree type.
      typedef int degree_type;
/// Alias for the numerical coefficient.
      typedef T numerical_type;
/// Alias for the rational coefficient.
      typedef mpq_class rational_type;
/// Alias for the container type.
      typedef std::valarray<expo_t> container_type;
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
          m.numerical_cf_=numerical_;
          m.rational_cf_=rational_;
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
          m.numerical_cf_=numerical_;
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
          m.rational_cf_=rational_;
        }
        const rational_type   &rational_;
      }
      ;
// Ctors.
/// Default constructor.
      explicit monomial_gmp_array():container_(0),numerical_cf_(0.),rational_cf_(0)
        {}
/// Constructor from size.
      explicit monomial_gmp_array(const size_t &n):container_(n),numerical_cf_(0.),rational_cf_(0)
        {}
/// Constructor from integer.
      explicit monomial_gmp_array(int n):container_(0),numerical_cf_(math::sgn(n)),
        rational_cf_((unsigned int)std::abs(n))
      {
        if (n==0)
        {
          numerical_cf_=numerical_type(0.);
          rational_cf_=1;
        }
      }
/// Constructor from numerical type.
      explicit monomial_gmp_array(const numerical_type &x):container_(0),
        numerical_cf_(x),rational_cf_(1)
        {}
/// Constructor from double.
      explicit monomial_gmp_array(const double &x):container_(0),
        numerical_cf_(x),rational_cf_(1)
        {}
      explicit monomial_gmp_array(const std::string &);
/// Constructor from numerical, rational and exponent container.
      explicit monomial_gmp_array(const numerical_type &num, const rational_type &rat, const container_type &cont)
        :container_(cont),numerical_cf_(num),rational_cf_(rat)
        {}
/// Copy constructor.
      monomial_gmp_array(const monomial_gmp_array &m):container_(m.container_),numerical_cf_(m.numerical_cf_),
        rational_cf_(m.rational_cf_)
        {}
/// Copy constructor from monomial with different numerical type.
      template <class U>
        monomial_gmp_array(const monomial_gmp_array<U> &m):container_(m.container()),
        numerical_cf_(m.numerical_cf()),rational_cf_(m.rational_cf())
        {}
/// Destructor.
      ~monomial_gmp_array()
        {}
// Getters
/// Get size of monomial (i.e., the number of variables).
      size_t width() const
      {
        return container_.size();
      }
/// Get const reference to the numerical coefficient.
      const numerical_type &numerical_cf() const
      {
        return numerical_cf_;
      }
/// Get reference to the numerical coefficient.
      numerical_type &numerical_cf()
      {
        return numerical_cf_;
      }
/// Get reference to the rational coefficient.
      const rational_type &rational_cf() const
      {
        return rational_cf_;
      }
      rational_type &rational_cf()
      {
        return rational_cf_;
      }
/// Get const reference to the container.
      const container_type &container() const
      {
        return container_;
      }
      container_type &container()
      {
        return container_;
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
      degree_type g_degree() const;
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
        p_assert(w==retval.width());
        const expo_t exponent=container_[n];
        if (exponent==0)
        {
// If the exponent is zero, partial derivative will be zero too.
          retval.numerical_cf()=numerical_type(0);
        }
        else
        {
// The rational cf of retval gets multiplied by the exponent.
          retval.rational_cf()=rational_cf_;
// Take abs because only numerical cf can be negative.
          retval.rational_cf()*=(expo_t)std::abs(exponent);
          if (exponent>0)
          {
            retval.numerical_cf()=numerical_cf_;
          }
          else
          {
            retval.numerical_cf()=-numerical_cf_;
          }
// Assign same exponents as this...
          retval.container_=container_;
// ... just modify the interested one.
          --retval.container_[n];
// TODO: place assert to make sure we don't go out expo_t range?
// This should be generalized, if we decide this way.
        }
      }
// Operators.
      template <class U>
        monomial_gmp_array &operator=(const monomial_gmp_array<U> &);
      bool operator<(const monomial_gmp_array &) const;
// End INTERFACE definition.
//-------------------------------------------------------
/// Reserve space for exponents.
/**
 * Postcondition: container is large enough to hold w elements. No assumptions can be made on
 * its contents.
 */
      void reserve(const size_t &w)
      {
        if (w!=width())
        {
          container_.resize(w);
        }
      }
    protected:
      template <class U>
        void basic_assignment(const monomial_gmp_array<U> &);
      template <class U>
        bool basic_comparison(const U &) const;
      bool is_symbolic() const;
    protected:
      container_type              container_;
      numerical_type              numerical_cf_;
      rational_type               rational_cf_;
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
    inline monomial_gmp_array<T>::monomial_gmp_array(const std::string &s):container_(0),numerical_cf_(0.),rational_cf_(0)
  {
    deque_string split_v;
    boost::split(split_v,s,boost::is_any_of(separator1_));
    if (split_v.size()==1)
    {
// In this case we try to cast to double.
      std::cout << "Trying to build monomial from numerical value." << std::endl;
      rational_cf_=1;
      numerical_cf_=utils::lexical_converter<numerical_type>(split_v[0]);
      return;
    }
    if (split_v.size()!=3)
    {
      std::cout << "Invalid number of elements in monomial_gmp_array string. Returning 0." << std::endl;
      return;
    }
    boost::trim(split_v[0]);
    numerical_cf_=utils::lexical_converter<numerical_type>(split_v[0]);
    boost::trim(split_v[1]);
    rational_cf_=utils::lexical_converter<rational_type>(split_v[1]);
    rational_cf_.canonicalize();
    if (rational_cf_<0)
    {
      rational_cf_=abs(rational_cf_);
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
    container_.resize(exponent_v.size());
    for (size_t i=0;i<exponent_v.size();++i)
    {
      if (exponent_v[i]!="")
      {
        container_[i]=utils::lexical_converter<expo_t>(exponent_v[i]);
      }
    }
  }

/// Print in plain format.
  template <class T>
    inline void monomial_gmp_array<T>::print_plain(std::ostream &out_stream, const vector_psym_p &) const
  {
    stream_manager::setup_print(out_stream);
    out_stream << numerical_cf_ << separator1_ << rational_cf_ << separator1_;
    out_stream << '[';
    for (size_t i=0;i<width();++i)
    {
      out_stream << container_[i];
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
    if (numerical_cf_.is_unity() && (is_symbolic() || rational_cf_ != 1))
    {
      if (numerical_cf_.is_negative())
      {
        tmp << std::string("-");
      }
    }
    else
    {
      tmp << numerical_cf_;
// Print cdot if the following fraction is really an integer - but not 1, because in that
// case it would not be printed.
      if (rational_cf_.get_den()==1 && rational_cf_.get_num()!=1)
      {
        tmp << std::string("\\cdot ");
      }
    }
    if (rational_cf_!=1)
    {
// If fraction is really an integer print only numerator.
      if (rational_cf_.get_den()==1)
      {
        tmp << rational_cf_.get_num();
      }
      else
      {
        tmp << "\\frac{" << rational_cf_.get_num() << "}{" << rational_cf_.get_den() << "}";
      }
    }
    for (size_t i=0;i<width();++i)
    {
      if (container_[i]==0)
        {}
        else if (container_[i]==1)
      {
        tmp << v[i]->name();
      }
      else
      {
// Trailing space needed to avoid problems with symbols whose name starts with a number.
        tmp << v[i]->name() << std::string("^{") << container_[i] << "} ";
      }
    }
    tmp << std::string("$");
// Print only if we have something to print.
    if (tmp.str()!=std::string("$$"))
    {
      out_stream << tmp.str();
    }
  }

/// Get degree of the monomial.
  template <class T>
    inline typename monomial_gmp_array<T>::degree_type monomial_gmp_array<T>::g_degree() const
  {
    degree_type retval=0;
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      retval+=container_[i];
    }
    return retval;
  }

/// Check whether monomial is symbolic or purely numerical.
  template <class T>
    inline bool monomial_gmp_array<T>::is_symbolic() const
  {
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      if (container_[i]!=0)
      {
        return true;
      }
    }
    return false;
  }

/// Hasher function.
  template <class T>
    inline size_t monomial_gmp_array<T>::hasher() const
  {
    size_t seed=0;
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      boost::hash_combine(seed,container_[i]);
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
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      if (container_[i]<m2.container_[i])
      {
        return true;
      }
      if (container_[i]>m2.container_[i])
      {
        return false;
      }
    }
    std::cout << "GWABBBBBBBB!!!!" << std::endl;
    std::exit(1);
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
      container_=m2.container();
      numerical_cf_=m2.numerical_cf();
      rational_cf_=m2.rational_cf();
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
    return (width()==w);
  }

/// Time evaluation.
  template <class T>
    inline typename monomial_gmp_array<T>::eval_type
    monomial_gmp_array<T>::t_eval(const double &t, const vector_psym_p &v) const
  {
    const size_t w=width();
    p_assert(v.size()==w);
    eval_type retval=numerical_cf_.value()*rational_cf_.get_d();
    for (size_t i=0;i<w;++i)
    {
// FIXME: use natural_pow or the like here, to speed up.
      retval*=std::pow(v[i]->t_eval(t),container_[i]);
    }
    return retval;
  }

/// Test for equality.
  template <class T>
    inline bool monomial_gmp_array<T>::operator==(const monomial_gmp_array &m2) const
  {
    return basic_comparison(m2);
  }

  template <class T>
    template <class U>
    inline bool monomial_gmp_array<T>::basic_comparison(const U &m2) const
  {
    const size_t w=width();
    p_assert(w==m2.width());
    for (size_t i=0;i<w;++i)
    {
      if (container_[i]!=m2.container_[i])
      {
        return false;
      }
    }
    return true;
  }

/// Increrase size so that final size is new_w.
  template <class T>
    inline void monomial_gmp_array<T>::increase_size(const size_t &new_w)
  {
    p_assert(new_w>=width());
    if (new_w>width())
    {
      container_type old=container_;
      container_.resize(new_w);
      const size_t old_w=old.size();
      for (size_t i=0;i<old_w;++i)
      {
        container_[i]=old[i];
      }
    }
  }

/// Prepend arguments.
  template <class T>
    inline void monomial_gmp_array<T>::prepend_args(const size_t &n)
  {
    if (n>0)
    {
      container_type old=container_;
      const size_t old_w=old.size();
      container_.resize(old_w+n);
      for (size_t i=0;i<old_w;++i)
      {
        container_[i+n]=old[i];
      }
    }
  }

/// Check whether monomial is zero.
// TODO: maybe this will have to include the handling of degree limit for polynomial variables.
  template <class T>
    inline bool monomial_gmp_array<T>::is_zero() const
  {
// We must do this way because if we get_d() immediately on rational_cf_ we may obtain zero for very
// small fractions, but the numerical_cf_ could counter-balance the smallness hence leading to a
// "non-zero" monomial.
    rational_type tmp(rational_cf_);
    tmp*=std::abs(numerical_cf_.value());
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
    out_m.numerical_cf_=numerical_cf_;
    out_m.numerical_cf_.mult_by_self(m2.numerical_cf());
    out_m.rational_cf_=rational_cf_;
    out_m.rational_cf_*=m2.rational_cf();
// Find big and small monomial_gmp_arrays.
    container_type const *big=&container_;
    container_type const *small=&m2.container();
    if (big->size()<small->size())
    {
      big=&m2.container();
      small=&container_;
    }
    out_m.container_.resize(big->size());
    const size_t max_w=big->size(), min_w=small->size();
    size_t i;
    for (i=0;i<min_w;++i)
    {
      out_m.container_[i]=(*big)[i]+(*small)[i];
    }
    for(;i<max_w;++i)
    {
      out_m.container_[i]=(*big)[i];
    }
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
  template <class T>
    inline monomial_gmp_array<T> monomial_gmp_array<T>::pow(int n_) const
  {
    p_assert(numerical_cf_.norm()>settings_manager::numerical_zero());
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
  }
}
#endif
