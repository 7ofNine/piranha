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

#ifndef PIRANHA_TRIG_ARRAY_H
#define PIRANHA_TRIG_ARRAY_H

#include "common_typedefs.h"    // For t_eval.
#include "trig_evaluator.h"

namespace piranha
{
/// Linear combination of arguments. Array version.
  class trig_array
  {
    public:
      size_t width() const;
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_array():private_flavour_(true)
        {}
/// Copy ctor.
      trig_array(const trig_array &ta):private_flavour_(ta.g_flavour()),private_container_(ta.private_container_)
        {}
      trig_array(const deque_string &);
      ~trig_array()
        {}
// Getters.
      mult_t multiplier(trig_size_t) const;
      size_t actual_width() const;
      bool &s_flavour()
      {
        return private_flavour_;
      }
      const bool &g_flavour() const
      {
        return private_flavour_;
      }
// I/O.
      void print_plain(std::ostream &, const vector_psym_p &) const;
      void print_latex(std::ostream &, const vector_psym_p &) const;
// Manip.
// FIXME: deprecate?
      void insert(trig_size_t,mult_t);
      void prepend_args(const size_t &);
      void append_args(const size_t &);
      void increase_size(const size_t &);
      void invert_sign();
// Probing.
      double freq(const vector_psym_p &) const;
      double phase(const vector_psym_p &) const;
      double t_eval(const double &, const vector_psym_p &) const;
      template <class TrigEvaluator>
        complex_double t_eval(TrigEvaluator &) const;
      short int sign() const;
      size_t hasher() const;
      bool is_zero() const;
      bool smaller(const size_t &) const;
      bool larger(const size_t &) const;
      bool compatible(const size_t &) const;
      size_t data_footprint() const;
      bool checkup(const size_t &) const;
      bool operator==(const trig_array &) const;
      bool operator<(const trig_array &) const;
// Math.
      void trigmult(const trig_array &, trig_array &, trig_array &) const;
      trig_array &operator=(const trig_array &);
      void operator*=(const mult_t &);
// End INTERFACE definition.
//-------------------------------------------------------
// Data members.
    private:
      bool                  private_flavour_;
      std::valarray<mult_t> private_container_;
  };

/// Ctor from piranha::deque_string.
  inline trig_array::trig_array(const deque_string &sd):private_flavour_(true),private_container_()
  {
    const size_t w=sd.size();
    if (w==0)
    {
      std::cout << "Warning: constructing empty trig_array." << std::endl;
      std::abort();
      return;
    }
// Now we know  w >= 1.
    private_container_.resize(w-1);
    for (size_t i=0;i<w-1;++i)
    {
      private_container_[i]=utils::lexical_converter<int>(sd[i]);
    }
// Take care of flavour.
    if (*sd.back().c_str()=='s')
    {
      s_flavour()=false;
    }
    else
    {
      s_flavour()=true;
    }
  }

// Getters implementations.
  inline mult_t trig_array::multiplier(trig_size_t index) const
  {
    return private_container_[index];
  }

  inline size_t trig_array::width() const
  {
    return private_container_.size();
  }

  inline size_t trig_array::actual_width() const
  {
    return width();
  }

// I/O implementations.
  inline void trig_array::print_plain(std::ostream &out_stream, const vector_psym_p &) const
  {
    stream_manager::setup_print(out_stream);
    for (size_t i=0;i<width();++i)
    {
      out_stream << private_container_[i] << stream_manager::data_separator();
    }
  }

  inline void trig_array::print_latex(std::ostream &out_stream, const vector_psym_p &v) const
  {
    stream_manager::setup_print(out_stream);
    bool first_one=true;
    std::string tmp("$");
    for (size_t i=0;i<width();++i)
    {
      if (private_container_[i]!=0)
      {
        if (private_container_[i]>0 && !first_one)
        {
          tmp.append("+");
        }
        if (private_container_[i]==-1)
        {
          tmp.append("-");
        }
        else if (private_container_[i]==1)
          {}
          else
        {
          tmp.append(boost::lexical_cast<std::string>(private_container_[i]));
        }
        tmp.append(v[i]->name());
        first_one=false;
      }
    }
    tmp.append("$");
// If we did not write anything erase math markers.
    if (tmp=="$$")
    {
      tmp.clear();
    }
    out_stream << tmp;
  }

// Manip implementations.
  inline void trig_array::insert(trig_size_t index, mult_t multiplier)
  {
    private_container_[index]=multiplier;
  }

  inline void trig_array::increase_size(const size_t &w)
  {
    p_assert(w>=width());
    if (w>width())
    {
      std::valarray<mult_t> old_private_container_(private_container_);
      private_container_.resize(w);
      const size_t old_w=old_private_container_.size();
      for (size_t i=0;i<old_w;++i)
      {
        private_container_[i]=old_private_container_[i];
      }
    }
  }

  inline void trig_array::append_args(const size_t &n)
  {
    increase_size(width()+n);
  }

  inline void trig_array::prepend_args(const size_t &n)
  {
    if (n>0)
    {
      std::valarray<mult_t> old_private_container_(private_container_);
      const size_t old_w=old_private_container_.size();
      private_container_.resize(old_w+n);
      for (size_t i=0;i<old_w;++i)
      {
        private_container_[i+n]=old_private_container_[i];
      }
    }
  }

// Probing implementations
/// Equality.
  inline bool trig_array::operator==(const trig_array &l2) const
  {
    const size_t w=width();
    p_assert(w==l2.width());
    for (size_t i=0;i<w;++i)
    {
      if (private_container_[i]!=l2.private_container_[i])
      {
        return false;
      }
    }
    return true;
  }

  inline bool trig_array::checkup(const size_t &w) const
  {
    if (w!=width())
    {
      std::cout << "Size mismatch in trig_array." << std::endl;
      return false;
    }
    return true;
  }

// FIXME: refactor these two to use a common low level function?
/// Frequency.
/**
 * Get the frequency of the linear combination, given a vector of piranha::psymbol pointers describing
 * the arguments.
 * @param[in] v vector of piranha::psymbol pointers.
 */
  inline double trig_array::freq(const vector_psym_p &v) const
  {
    p_assert(v.size()==width());
    double retval=0.;
    for (size_t i=0;i<width();++i)
    {
// We must be sure that there actually is a freq in every symbol we are going to use.
      if (v[i]->poly_eval().size()>1)
      {
        retval+=private_container_[i]*v[i]->poly_eval()[1];
      }
    }
    return retval;
  }

/// Phase.
/**
 * Get the phase of the linear combination, given a vector of piranha::psymbol pointers describing the
 * arguments.
 * @param[in] v vector of piranha::psymbol pointers.
 */
  inline double trig_array::phase(const vector_psym_p &v) const
  {
    p_assert(v.size()==width());
    double retval=0.;
    for (size_t i=0;i<width();++i)
    {
// We must be sure that there actually is a phase in every symbol we are going to use.
      if (v[i]->poly_eval().size()>0)
      {
        retval+=private_container_[i]*v[i]->poly_eval()[0];
      }
    }
    return retval;
  }

/// Time evaluation of arguments.
/**
 * Returns the value assumed by the linear combination of arguments at time t.
 * @param[in] t double time of the evaluation.
 * @param[in] v vector of piranha::psymbol pointers.
 */
  inline double trig_array::t_eval(const double &t, const vector_psym_p &v) const
  {
    p_assert(v.size()==width());
    double retval=0.;
    for (size_t i=0;i<v.size();++i)
    {
      if (private_container_[i]!=0)
      {
        retval+=private_container_[i]*v[i]->t_eval(t);
      }
    }
    return retval;
  }

/// Time evaluation of complex exponential of the arguments.
/**
 * Returns the complex exponential of the linear combination of arguments at time t.
 * Uses a TrigEvaluator objet which contains a cache of the complex exponentials of arguments.
 * @param[in] te piranha::trig_evaluator containing a cache of complex exponentials of arguments.
 */
  template <class TrigEvaluator>
    inline complex_double trig_array::t_eval(TrigEvaluator &te) const
  {
    p_assert(te.width()==width());
    const size_t w=width();
    complex_double retval(1.);
    for (size_t i=0;i<w;++i)
    {
      if (private_container_[i]!=0)
      {
        retval*=te.request_complexp(i,private_container_[i]);
      }
    }
    return retval;
  }

/// Sign.
/**
 * Return the sign of the first non-zero element of the combination. Zero is considered positive.
 */
  inline short int trig_array::sign() const
  {
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      if (private_container_[i]>0)
      {
        return 1;
      }
      if (private_container_[i]<0)
      {
        return -1;
      }
    }
    return 1;
  }

/// Data footprint.
/**
 * Returns the memory occupied by the data members.
 */
  inline size_t trig_array::data_footprint() const
  {
    return (width()*sizeof(mult_t));
  }

/// Operator less than.
/**
 * Needed, for example, in norm-based inidces for piranha::base_pseries, where there could
 * be terms with equal norms. If that happens, this operator is used to order terms.
 */
  inline bool trig_array::operator<(const trig_array &l2) const
  {
    size_t w=width();
    p_assert(w==l2.width());
    for (size_t i=0;i<w;++i)
    {
      if (private_container_[i]<l2.private_container_[i])
      {
        return true;
      }
      else if (private_container_[i]>l2.private_container_[i])
      {
        return false;
      }
    }
    return false;
  }

  inline size_t trig_array::hasher() const
  {
    size_t seed=0;
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      boost::hash_combine(seed,private_container_[i]);
    }
    return seed;
  }

  inline bool trig_array::is_zero() const
  {
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      if (private_container_[i]!=0)
      {
        return false;
      }
    }
    return true;
  }

  inline bool trig_array::smaller(const size_t &n) const
  {
    return (width()<n);
  }

  inline bool trig_array::larger(const size_t &n) const
  {
    return (width()>n);
  }

  inline bool trig_array::compatible(const size_t &n) const
  {
    return (width()==n);
  }

/// Multiply by a piranha::mult_t.
  inline void trig_array::operator*=(const mult_t &n)
  {
    const size_t w=width();
    for (size_t i=0;i<w;++i)
    {
      private_container_[i]*=n;
    }
  }

/// Invert sign.
  inline void trig_array::invert_sign()
  {
    *this*=-1;
  }

/// Assignment operator.
/**
 * After the assignment the data members of the two classes must be equal, both from a mathematical
 * point of view and with respect to the computer representation. Assignment to a larger trig_array
 * is allowed, assignment to smaller results in assertion failure.
 * @param[in] l2 right-hand side piranha::trig_array.
 */
  inline trig_array &trig_array::operator=(const trig_array &l2)
  {
    if (&l2!=this)
    {
      increase_size(l2.width());
      s_flavour()=l2.g_flavour();
      private_container_=l2.private_container_;
    }
    return *this;
  }

/// Multiplication.
/**
 * Multiplication of two trigonometric functions using Werner's formulas, i.e.
 * \f[
 * C\cos\alpha\cdot\cos\beta=
 * \frac{C}{2} \cos \left( \alpha - \beta \right) + \frac{C}{2} \cos \left( \alpha + \beta \right)
 * \f]
 * and the likes. Notice that in the first return value always goes the \f$ \alpha - \beta \f$ term
 * and in the second one always goes \f$ \alpha + \beta \f$ one.
 * Please also note that no assumptions are made with respect to return values' content (e.g., it is not guaranteed
 * that return values are empty).
 * @param[in] l2 factor.
 * @param[out] ret1 first return value.
 * @param[out] ret2 second return value.
 */
// NOTE: we are not using here a general version of vector addition/subtraction
// because this way we can do two operations (+ and -) every cycle. This is a performance
// critical part, so the optimization should be worth the hassle.
  inline void trig_array::trigmult(const trig_array &l2, trig_array &ret1, trig_array &ret2) const
  {
    const size_t max_w=width(), min_w=l2.width();
// Assert widths, *this should always come from a regular Poisson series, and its width should hence be
// already adjusted my merge_args in multiplication routines.
    p_assert(max_w >= min_w);
    size_t i;
    ret1.increase_size(max_w);
    ret2.increase_size(max_w);
    for (i=0;i<min_w;++i)
    {
      ret1.private_container_[i]=private_container_[i]-l2.private_container_[i];
      ret2.private_container_[i]=private_container_[i]+l2.private_container_[i];
    }
    for (;i<max_w;++i)
    {
      ret1.private_container_[i]=private_container_[i];
      ret2.private_container_[i]=private_container_[i];
    }
  }

/// Overload of hash_value function for piranha::trig_array.
/**
 * To be used in piranha::base_pseries for the hashed index.
 */
  inline size_t hash_value(const trig_array &ta)
  {
    return ta.hasher();
  }
}
#endif
