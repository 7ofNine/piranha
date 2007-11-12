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

#ifndef PIRANHA_BASE_TRIG_ARRAY_H
#define PIRANHA_BASE_TRIG_ARRAY_H

#include <boost/integer.hpp>
#include <boost/static_assert.hpp>

#include "../bits/common_typedefs.h"    // For t_eval.
#include "../bits/trig_evaluator.h"

namespace piranha
{
/// Base class for dense trigonometric array.
/**
 * Inherited by specialized trigonometric array classes like piranha::trig_array and trig_fixed_array.
 */
  template <int Bits, class Derived>
    class base_trig_array
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    public:
      typedef typename boost::int_t<Bits>::fast value_type;
// Ctors.
/// Default ctor.
      base_trig_array():private_flavour_(true)
        {}
/// Copy ctor.
      base_trig_array(const base_trig_array &t):private_flavour_(t.g_flavour())
        {}
      ~base_trig_array()
        {}
// Getters.
      value_type at(const trig_size_t &n) const
      {
        if (static_cast<const Derived *>(this)->g_width() <= n)
        {
          std::cout << "Warning, trying to extract out-of-boundaries element from trigonometric array." << std::endl;
// TODO: drop this abort later when we log properly.
          std::abort();
          return 0;
        }
        return static_cast<const Derived *>(this)->g_container()[n];
      }
      size_t actual_width() const
      {
        return static_cast<const Derived *>(this)->g_width();
      }
      bool &s_flavour()
      {
        return private_flavour_;
      }
      const bool &g_flavour() const
      {
        return private_flavour_;
      }
// I/O.
      void print_plain(std::ostream &out_stream, const vector_psym_p &) const
      {
        stream_manager::setup_print(out_stream);
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          out_stream << (int)static_cast<const Derived *>(this)->g_container()[i] << stream_manager::data_separator();
        }
      }
      void print_latex(std::ostream &out_stream, const vector_psym_p &v) const
      {
        stream_manager::setup_print(out_stream);
        bool first_one=true;
        std::string tmp("$");
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i]!=0)
          {
            if (static_cast<const Derived *>(this)->g_container()[i]>0 && !first_one)
            {
              tmp.append("+");
            }
            if (static_cast<const Derived *>(this)->g_container()[i]==-1)
            {
              tmp.append("-");
            }
            else if (static_cast<const Derived *>(this)->g_container()[i]==1)
              {}
              else
            {
              tmp.append(boost::lexical_cast<std::string>((int)static_cast<const Derived *>(this)->g_container()[i]));
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
      void invert_sign()
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          static_cast<const Derived *>(this)->s_container()[i]=
          -static_cast<const Derived *>(this)->g_container()[i];
        }
      }
/// Assign vector of multipliers.
      template <class T>
        void assign_mult_vector(const T &v)
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        p_assert(v.size() == w);
        for (trig_size_t i=0;i<w;++i)
        {
          static_cast<const Derived *>(this)->s_container()[i]=v[i];
        }
      }
      template <class DerivedPs>
        double density(const DerivedPs &p) const
      {
        size_t tmp=0;
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i] != 0)
          {
            ++tmp;
          }
        }
        return ((double)tmp)/w;
      }
// TODO: refactor these two to use a common low level function?
// TODO: is there some caching mechanism that can be used here?
/// Frequency.
/**
 * Get the frequency of the linear combination, given a vector of piranha::psymbol pointers describing
 * the arguments.
 * @param[in] v vector of piranha::psymbol pointers.
 */
      double freq(const vector_psym_p &v) const
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        p_assert(v.size() == w);
        double retval=0.;
        for (trig_size_t i=0;i<w;++i)
        {
// We must be sure that there actually is a freq in every symbol we are going to use.
          if (v[i]->poly_eval().size()>1)
          {
            retval+=static_cast<const Derived *>(this)->g_container()[i]*v[i]->poly_eval()[1];
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
      double phase(const vector_psym_p &v) const
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        p_assert(v.size()==w);
        double retval=0.;
        for (trig_size_t i=0;i<w;++i)
        {
// We must be sure that there actually is a phase in every symbol we are going to use.
          if (v[i]->poly_eval().size()>0)
          {
            retval+=static_cast<const Derived *>(this)->g_container()[i]*v[i]->poly_eval()[0];
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
      double t_eval(const double &t, const vector_psym_p &v) const
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        p_assert(v.size()==w);
        double retval=0.;
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i]!=0)
          {
            retval+=static_cast<const Derived *>(this)->g_container()[i]*v[i]->t_eval(t);
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
        complex_double t_eval(TrigEvaluator &te) const
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        p_assert(te.width() == w);
        complex_double retval(1.);
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i]!=0)
          {
            retval*=te.request_complexp(i,static_cast<const Derived *>(this)->g_container()[i]);
          }
        }
        return retval;
      }
/// Sign.
/**
 * Return the sign of the first non-zero element of the combination. Zero is considered positive.
 */
      short int sign() const
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i]>0)
          {
            return 1;
          }
          if (static_cast<const Derived *>(this)->g_container()[i]<0)
          {
            return -1;
          }
        }
        return 1;
      }
      size_t hasher() const
      {
        size_t seed=g_flavour();
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (size_t i=0;i<w;++i)
        {
          boost::hash_combine(seed,static_cast<const Derived *>(this)->g_container()[i]);
        }
        return seed;
      }
      bool is_zero() const
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i]!=0)
          {
            return false;
          }
        }
        return true;
      }
      bool smaller(const size_t &n) const
      {
        return (static_cast<const Derived *>(this)->g_width() < n);
      }
      bool larger(const size_t &n) const
      {
        return (static_cast<const Derived *>(this)->g_width() > n);
      }
      bool compatible(const size_t &n) const
      {
        return (static_cast<const Derived *>(this)->g_width() == n);
      }
    protected:
      void assignment_operator(const Derived &t2)
      {
        if (static_cast<const Derived *>(this) != &t2)
        {
          private_flavour_=t2.private_flavour_;
          static_cast<const Derived *>(this)->assignment(t2);
        }
      }
      bool equality_test(const Derived &t2) const
      {
        p_assert(static_cast<const Derived *>(this)->g_width() == t2.g_width());
        if (g_flavour() != t2.g_flavour())
        {
          return false;
        }
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i] != t2.g_container()[i])
          {
            return false;
          }
        }
        return true;
      }
/// Less than.
      bool less_than(const Derived &t2) const
      {
        if (g_flavour() < t2.g_flavour())
        {
          return true;
        }
        else if (g_flavour() > t2.g_flavour())
        {
          return false;
        }
// TODO: maybe when unrolling here we have to substitute direct "return" calls with retval and break statements.
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          if (static_cast<const Derived *>(this)->g_container()[i] < t2.g_container()[i])
          {
            return true;
          }
          else if (static_cast<const Derived *>(this)->g_container()[i] > t2.g_container()[i])
          {
            return false;
          }
        }
        return false;
      }
/// Multiply by a int.
      void mult_by_int(const int &n)
      {
        const trig_size_t w=static_cast<const Derived *>(this)->g_width();
        for (trig_size_t i=0;i<w;++i)
        {
          static_cast<const Derived *>(this)->s_container()[i]*=n;
        }
      }
    private:
      bool  private_flavour_;
  };

/// Overload of hash_value function for piranha::base_trig_array.
/**
 * To be used in piranha::base_pseries for the hashed index.
 */
  template <int Bits, class Derived>
    inline size_t hash_value(const base_trig_array<Bits,Derived> &t)
  {
    return static_cast<const Derived *>(&t)->hasher();
  }
}

#endif