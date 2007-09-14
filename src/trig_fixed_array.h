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

#ifndef PIRANHA_TRIG_FIXED_ARRAY_H
#define PIRANHA_TRIG_FIXED_ARRAY_H

#include <boost/static_assert.hpp>
#include <boost/algorithm/minmax.hpp>
#include <cstring>

#include "common_typedefs.h"    // For t_eval.
#include "trig_evaluator.h"

namespace piranha
{
/// Linear combination of arguments. Fixed size array version.
  template <int Dim>
    class trig_fixed_array
  {
// Check that dimension is sane.
      BOOST_STATIC_ASSERT(Dim > 0);
      BOOST_STATIC_ASSERT(Dim < 100);
    public:
// Start INTERFACE definition.
//-------------------------------------------------------
// Ctors.
/// Default ctor.
      trig_fixed_array():private_flavour_(true)
        {}
/// Copy ctor.
      trig_fixed_array(const trig_fixed_array &ta):private_flavour_(ta.g_flavour())
        {
          memcpy((void *)private_container_,(const void *)ta.private_container_,sizeof(mult_t)*Dim);
        }
      trig_fixed_array(const deque_string &);
      ~trig_fixed_array()
        {}
// Getters.
      mult_t at(trig_size_t) const;
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
      void prepend_args(const size_t &);
      void append_args(const size_t &);
      void increase_size(const size_t &);
      void invert_sign();
      template <class T>
        void assign_mult_vector(const T &);
// Probing.
      template <class DerivedPs>
        double density(const DerivedPs &p) const
      {
// TODO: unroll.
        size_t tmp=0;
        for (size_t i=0;i<Dim;++i)
        {
          if (private_container_[i] != 0)
          {
            ++tmp;
          }
        }
        return ((double)tmp)/Dim;
      }
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
      bool operator==(const trig_fixed_array &) const;
      bool operator<(const trig_fixed_array &) const;
// Math.
      void trigmult(const trig_fixed_array &, trig_fixed_array &, trig_fixed_array &) const;
      trig_fixed_array &operator=(const trig_fixed_array &);
      void operator*=(const mult_t &);
// End INTERFACE definition.
//-------------------------------------------------------
    private:
// Data members.
    private:
      bool    private_flavour_;
      mult_t  private_container_[Dim];
  };

/// Ctor from piranha::deque_string.
  template <int Dim>
    inline trig_fixed_array<Dim>::trig_fixed_array(const deque_string &sd):private_flavour_(true)
  {
    const size_t w=sd.size();
    if (w == 0)
    {
      std::cout << "Warning: not enough elements to construct trig_fixed_array." << std::endl;
      std::abort();
      return;
    }
// Now w >= 1.
    if ((w-1) != Dim)
    {
      std::cout << "Warning: wrong size for trig_fixed_array in ctor from string." << std::endl;
// TODO: Here we continue really, just remains as debug.
      std::abort();
    }
    size_t i;
    for (i=0;i<boost::minmax((size_t)Dim,w-1).get<0>();++i)
    {
      private_container_[i]=utils::lexical_converter<int>(sd[i]);
    }
    for (;i<Dim;++i)
    {
      private_container_[i]=0;
    }
// Take care of flavour.
    if (*sd.back().c_str()=='s')
    {
      s_flavour()=false;
    }
  }

// Getters implementations.
  template <int Dim>
    inline mult_t trig_fixed_array<Dim>::at(trig_size_t index) const
  {
    return private_container_[index];
  }

  template <int Dim>
    inline size_t trig_fixed_array<Dim>::actual_width() const
  {
    return Dim;
  }

// I/O implementations.
// TODO: place asserts against vector_psym_p width? Also in other array trigs?
  template <int Dim>
    inline void trig_fixed_array<Dim>::print_plain(std::ostream &out_stream, const vector_psym_p &) const
  {
    stream_manager::setup_print(out_stream);
    for (size_t i=0;i<Dim;++i)
    {
      out_stream << private_container_[i] << stream_manager::data_separator();
    }
  }

  template <int Dim>
    inline void trig_fixed_array<Dim>::print_latex(std::ostream &out_stream, const vector_psym_p &v) const
  {
    stream_manager::setup_print(out_stream);
    bool first_one=true;
    std::string tmp("$");
    for (size_t i=0;i<Dim;++i)
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
/// Assign vector of multipliers.
// TODO: unroll? memcpy?
  template <int Dim>
    template <class T>
    inline void trig_fixed_array<Dim>::assign_mult_vector(const T &v)
  {
    p_assert(v.size() == Dim);
    for (size_t i=0;i<Dim;++i)
    {
      private_container_[i]=v[i];
    }
  }

  template <int Dim>
    inline void trig_fixed_array<Dim>::increase_size(const size_t &w)
  {
    p_assert(w==Dim);
  }

  template <int Dim>
    inline void trig_fixed_array<Dim>::append_args(const size_t &n)
  {
    p_assert(false);
  }

  template <int Dim>
    inline void trig_fixed_array<Dim>::prepend_args(const size_t &n)
  {
    p_assert(false);
  }

// Probing implementations
/// Equality.
// TODO: unroll.
  template <int Dim>
    inline bool trig_fixed_array<Dim>::operator==(const trig_fixed_array &l2) const
  {
    if (g_flavour() != l2.g_flavour())
    {
      return false;
    }
    for (size_t i=0;i<Dim;++i)
    {
      if (private_container_[i] != l2.private_container_[i])
      {
        return false;
      }
    }
    return true;
  }

  template <int Dim>
    inline bool trig_fixed_array<Dim>::checkup(const size_t &w) const
  {
    return true;
  }

// FIXME: refactor these two to use a common low level function?
// TODO: is there some caching mechanism that can be used here?
/// Frequency.
/**
 * Get the frequency of the linear combination, given a vector of piranha::psymbol pointers describing
 * the arguments.
 * @param[in] v vector of piranha::psymbol pointers.
 */
  template <int Dim>
    inline double trig_fixed_array<Dim>::freq(const vector_psym_p &v) const
  {
    p_assert(v.size()==Dim);
    double retval=0.;
    for (size_t i=0;i<Dim;++i)
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
  template <int Dim>
    inline double trig_fixed_array<Dim>::phase(const vector_psym_p &v) const
  {
    p_assert(v.size()==Dim);
    double retval=0.;
    for (size_t i=0;i<Dim;++i)
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
// TODO: unroll.
  template <int Dim>
    inline double trig_fixed_array<Dim>::t_eval(const double &t, const vector_psym_p &v) const
  {
    p_assert(v.size()==Dim);
    double retval=0.;
    for (size_t i=0;i<Dim;++i)
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
// TODO: unroll.
  template <int Dim>
    template <class TrigEvaluator>
    inline complex_double trig_fixed_array<Dim>::t_eval(TrigEvaluator &te) const
  {
    p_assert(te.width()==Dim);
    complex_double retval(1.);
    for (size_t i=0;i<Dim;++i)
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
// TODO: unroll.
  template <int Dim>
    inline short int trig_fixed_array<Dim>::sign() const
  {
    for (size_t i=0;i<Dim;++i)
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
  template <int Dim>
    inline size_t trig_fixed_array<Dim>::data_footprint() const
  {
    return sizeof(trig_fixed_array);
  }

/// Operator less than.
/**
 * Needed, for example, in norm-based inidces for piranha::base_pseries, where there could
 * be terms with equal norms. If that happens, this operator is used to order terms.
 */
// TODO: unroll.
  template <int Dim>
    inline bool trig_fixed_array<Dim>::operator<(const trig_fixed_array &l2) const
  {
    if (g_flavour() < l2.g_flavour())
    {
      return true;
    }
    else if (g_flavour() > l2.g_flavour())
    {
      return false;
    }
// TODO: maybe when unrolling here we have to substitute direct "return" calls with retval and break statements.
    for (size_t i=0;i<Dim;++i)
    {
      if (private_container_[i] < l2.private_container_[i])
      {
        return true;
      }
      else if (private_container_[i] > l2.private_container_[i])
      {
        return false;
      }
    }
    return false;
  }

// TODO: unroll.
  template <int Dim>
    inline size_t trig_fixed_array<Dim>::hasher() const
  {
    size_t seed=g_flavour();
    for (size_t i=0;i<Dim;++i)
    {
      boost::hash_combine(seed,private_container_[i]);
    }
    return seed;
  }

// TODO: unroll.
  template <int Dim>
    inline bool trig_fixed_array<Dim>::is_zero() const
  {
    for (size_t i=0;i<Dim;++i)
    {
      if (private_container_[i]!=0)
      {
        return false;
      }
    }
    return true;
  }

  template <int Dim>
    inline bool trig_fixed_array<Dim>::smaller(const size_t &n) const
  {
    return (Dim<n);
  }

  template <int Dim>
    inline bool trig_fixed_array<Dim>::larger(const size_t &n) const
  {
    return (Dim>n);
  }

  template <int Dim>
    inline bool trig_fixed_array<Dim>::compatible(const size_t &n) const
  {
    return (Dim==n);
  }

/// Multiply by a piranha::mult_t.
// TODO: unroll.
  template <int Dim>
    inline void trig_fixed_array<Dim>::operator*=(const mult_t &n)
  {
    for (size_t i=0;i<Dim;++i)
    {
      private_container_[i]*=n;
    }
  }

/// Invert sign.
// TODO: unroll.
  template <int Dim>
    inline void trig_fixed_array<Dim>::invert_sign()
  {
    for (size_t i=0;i<Dim;++i)
    {
      private_container_[i]=-private_container_[i];
    }
  }

/// Assignment operator.
/**
 * After the assignment the data members of the two classes must be equal, both from a mathematical
 * point of view and with respect to the computer representation. Assignment to a larger trig_fixed_array
 * is allowed, assignment to smaller results in assertion failure.
 * @param[in] l2 right-hand side piranha::trig_fixed_array.
 */
// TODO: unroll or memcpy?
  template <int Dim>
    inline trig_fixed_array<Dim> &trig_fixed_array<Dim>::operator=(const trig_fixed_array &l2)
  {
    if (&l2!=this)
    {
      s_flavour()=l2.g_flavour();
      memcpy((void *)private_container_,(const void *)l2.private_container_,sizeof(mult_t)*Dim);
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
  template <int Dim>
    inline void trig_fixed_array<Dim>::trigmult(const trig_fixed_array &l2, trig_fixed_array &ret1, trig_fixed_array &ret2) const
  {
    for (size_t i=0;i<Dim;++i)
    {
      ret1.private_container_[i]=private_container_[i]-l2.private_container_[i];
      ret2.private_container_[i]=private_container_[i]+l2.private_container_[i];
    }
  }

/// Overload of hash_value function for piranha::trig_fixed_array.
/**
 * To be used in piranha::base_pseries for the hashed index.
 */
  template <int Dim>
    inline size_t hash_value(const trig_fixed_array<Dim> &ta)
  {
    return ta.hasher();
  }
}
#endif
