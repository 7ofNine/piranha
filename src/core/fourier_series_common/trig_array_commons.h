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

#ifndef PIRANHA_TRIG_ARRAY_COMMONS_H
#define PIRANHA_TRIG_ARRAY_COMMONS_H

#include <boost/integer.hpp>
#include <boost/static_assert.hpp>
#include <complex>
#include <string>
#include <vector>

#include "../common_typedefs.h" // For t_eval, max_fast_int and layout.
#include "../psymbol.h"
#include "trig_evaluator.h"

#define derived_const_cast (static_cast<Derived const *>(this))
#define derived_cast (static_cast<Derived *>(this))

namespace piranha
{
  /// Common class for dense trigonometric array.
  /**
   * Intended to add specific methods to plain arrays for the manipulation of trigonometric
   * parts in Poisson series.
   */
  template <class Derived>
    class trig_array_commons
  {
    public:
      // I/O.
      template <class ArgsTuple>
        void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const
      {
        // We assert like this because we want to make sure we don't go out of boundaries,
        // and because in case of fixed-width we may have smaller size of v wrt to "real" size.
        p_assert(args_tuple.template get<Derived::position>().size() <= derived_const_cast->size());
        derived_const_cast->print(out_stream);
        out_stream << derived_const_cast->separator;
        switch (derived_const_cast->flavour())
        {
          case true:
            out_stream << "c";
            break;
          case false:
            out_stream << "s";
        }
      }
      void print_latex(std::ostream &out_stream, const vector_psym_p &v) const
      {
        const size_t w=v.size();
        p_assert(w <= derived_const_cast->size())
          stream_manager::setup_print(out_stream);
        switch (derived_const_cast->flavour())
        {
          case true:
            out_stream << "c&";
            break;
          case false:
            out_stream << "s&";
        }
        bool first_one=true;
        std::string tmp("$");
        for (size_t i=0;i < w;++i)
        {
          if ((*derived_const_cast)[i] != 0)
          {
            if ((*derived_const_cast)[i] > 0 and !first_one)
            {
              tmp.append("+");
            }
            if ((*derived_const_cast)[i] == -1)
            {
              tmp.append("-");
            }
            else if ((*derived_const_cast)[i] == 1)
              {}
              else
            {
              tmp.append(boost::lexical_cast<std::string>((max_fast_int)(*derived_const_cast)[i]));
            }
            tmp.append(v[i]->name());
            first_one=false;
          }
        }
        tmp.append("$");
        // If we did not write anything erase math markers.
        if (tmp == "$$")
        {
          tmp.clear();
        }
        out_stream << tmp;
      }
      void invert_sign()
      {
        const size_t w=derived_const_cast->size();
        for (size_t i=0;i < w;++i)
        {
          (*derived_cast)[i]*=-1;
        }
      }
      /// Frequency.
      /**
       * Get the frequency of the linear combination, given a vector of piranha::psymbol pointers describing
       * the arguments.
       * @param[in] v vector of piranha::psymbol pointers.
       */
      template <class ArgsTuple>
        double freq(const ArgsTuple &args_tuple) const {return combined_time_eval<1>(args_tuple);}
      /// Phase.
      /**
       * Get the phase of the linear combination, given a vector of piranha::psymbol pointers describing the
       * arguments.
       * @param[in] v vector of piranha::psymbol pointers.
       */
      template <class ArgsTuple>
        double phase(const ArgsTuple &args_tuple) const {return combined_time_eval<0>(args_tuple);}
      /// Time evaluation of arguments.
      /**
       * Returns the value assumed by the linear combination of arguments at time t.
       * @param[in] t double time of the evaluation.
       * @param[in] v vector of piranha::psymbol pointers.
       */
      template <class ArgsTuple>
        double t_eval(const double &t, const ArgsTuple &args_tuple) const
      {
        const size_t w=args_tuple.template get<Derived::position>().size();
        p_assert(w <= derived_const_cast->size());
        double retval=0.;
        for (size_t i=0;i < w;++i)
        {
          if ((*derived_const_cast)[i] != 0)
          {
            retval+=(*derived_const_cast)[i]*args_tuple.template get<Derived::position>()[i]->t_eval(t);
          }
        }
        switch (derived_const_cast->flavour())
        {
          case true:
            return std::cos(retval);
          default:
            return std::sin(retval);
        }
      }
      /// Time evaluation of complex exponential of the arguments.
      /**
       * Returns the real or imaginary part (depending on flavour) of the complex exponential of the
       * linear combination of arguments at time t.
       * Uses a piranha::trig_evaluator object which contains a cache of the complex exponentials of arguments.
       * @param[in] te piranha::trig_evaluator containing a cache of complex exponentials of arguments.
       */
      template <class TrigEvaluator>
        double t_eval(TrigEvaluator &te) const
      {
        const size_t w=te.width();
        p_assert(w <= derived_const_cast->size());
        std::complex<double> retval(1.);
        for (size_t i=0;i < w;++i)
        {
          if ((*derived_const_cast)[i] != 0)
          {
            retval*=te.request_complexp(i,(*derived_const_cast)[i]);
          }
        }
        switch (derived_const_cast->flavour())
        {
          case true:
            return retval.real();
          default:
            return retval.imag();
        }
      }
      /// Sign.
      /**
       * Return the sign of the first non-zero element of the combination. Zero is considered positive.
       */
      short int sign() const
      {
        const size_t w=derived_const_cast->size();
        for (size_t i=0;i < w;++i)
        {
          // TODO: use switch?
          if ((*derived_const_cast)[i] > 0)
          {
            return 1;
          }
          if ((*derived_const_cast)[i] < 0)
          {
            return -1;
          }
        }
        return 1;
      }
      // All multipliers are zero and flavour is sine.
      template <class ArgsTuple>
        bool is_ignorable(const ArgsTuple &) const
      {
        return (derived_const_cast->is_zero() and !derived_const_cast->flavour());
      }
      /// Equality test.
      bool operator==(const Derived &t2) const
      {
        return (derived_const_cast->flavour() == t2.flavour() and derived_const_cast->equal_to(t2));
      }
      /// Less than.
      bool operator<(const Derived &t2) const
      {
        if (derived_const_cast->flavour() < t2.flavour())
        {
          return true;
        }
        else if (derived_const_cast->flavour() > t2.flavour())
        {
          return false;
        }
        const size_t w=derived_const_cast->size();
        for (size_t i=0;i<w;++i)
        {
          if ((*derived_const_cast)[i] < t2[i])
          {
            return true;
          }
          else if ((*derived_const_cast)[i] > t2[i])
          {
            return false;
          }
        }
        return false;
      }
      size_t hash_value() const
      {
        size_t retval = derived_const_cast->hasher();
        boost::hash_combine(retval,derived_const_cast->flavour());
        return retval;
      }
      /// Multiply by an integer.
      void mult_by_int(const int &n)
      {
        const size_t w=derived_const_cast->size();
        for (size_t i=0;i < w;++i)
        {
          (*derived_cast)[i]*=n;
        }
      }
    protected:
      trig_array_commons() {}
      trig_array_commons(const std::string &s)
      {
        typedef typename Derived::value_type value_type;
        std::vector<std::string> sd;
        boost::split(sd,s,boost::is_any_of(std::string(1,derived_const_cast->separator)));
        // TODO: check here that we are not loading too many multipliers, outside trig_size_t range.
        // TODO: do it everywhere!
        const size_t w=sd.size();
        if (w == 0)
        {
          std::cout << "Warning: constructing empty trig_array." << std::endl;
          std::abort();
          return;
        }
        // Now we know  w >= 1.
        derived_cast->resize(w-1);
        for (size_t i=0;i < w-1;++i)
        {
          (*derived_cast)[i]=utils::lexical_converter<value_type>(sd[i]);
        }
        // Take care of flavour.
        if (*sd.back().c_str() == 's')
        {
          derived_cast->flavour()=false;
        }
        else
        {
          derived_cast->flavour()=true;
        }
      }
    private:
      // NOTICE: is there some caching mechanism that can be used here?
      template <int N, class ArgsTuple>
        double combined_time_eval(const ArgsTuple &args_tuple) const
      {
        BOOST_STATIC_ASSERT(N >= 0);
        const size_t w=args_tuple.template get<Derived::position>().size();
        p_assert(w <= derived_const_cast->size());
        double retval=0.;
        for (size_t i=0;i < w;++i)
        {
          // We must be sure that there actually is component N in every symbol we are going to use.
          if (args_tuple.template get<Derived::position>()[i]->time_eval().size() > N)
          {
            retval+=(*derived_const_cast)[i]*args_tuple.template get<Derived::position>()[i]->time_eval()[N];
          }
        }
        return retval;
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
