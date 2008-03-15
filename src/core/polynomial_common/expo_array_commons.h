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

#ifndef PIRANHA_EXPO_ARRAY_COMMONS_H
#define PIRANHA_EXPO_ARRAY_COMMONS_H

#include <boost/algorithm/string/split.hpp>
#include <boost/static_assert.hpp>
#include <cmath> // For std::pow, most likely temporary.
#include <iostream>
#include <string>
#include <vector>

#define derived_const_cast (static_cast<Derived const *>(this))
#define derived_cast (static_cast<Derived *>(this))

namespace piranha
{
  /// Common class for dense exponent array.
  /**
   * Intended to extend piranha::int_array for the manipulation of exponent
   * parts in polynomials.
   */
  template <class Derived>
    class expo_array_commons
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
      }
      void print_latex(std::ostream &out_stream, const vector_psym_p &v) const
      {
        // TODO: implement.
      }
      template <class ArgsTuple>
        double t_eval(const double &t, const ArgsTuple &args_tuple) const
      {
        const size_t w=args_tuple.template get<Derived::position>().size();
        p_assert(w <= derived_const_cast->size());
        double retval=1.;
        for (size_t i=0;i < w;++i)
        {
          retval*=std::pow(args_tuple.template get<Derived::position>()[i]->t_eval(t),(*derived_const_cast)[i]);
        }
        return retval;
      }
      /// Am I ignorabled?
      /**
       * A product of integer powers of literal variables is never zero, hence never ignorable.
       */
      template <class ArgsTuple>
        bool is_ignorable(const ArgsTuple &) const
      {
        return false;
      }
      /// Equality test.
      bool operator==(const Derived &e2) const
      {
        return derived_const_cast->equal_to(e2);
      }
      /// Calculate hash value.
      size_t hash_value() const
      {
        return derived_const_cast->hasher();
      }
      /// Return the total degree of the exponents array.
      int get_degree() const
      {
        int retval = 0;
        for (typename Derived::size_type i=0; i < derived_const_cast->m_size; ++i)
        {
          retval+=(*derived_const_cast)[i];
        }
        return retval;
      }
    protected:
      expo_array_commons() {}
      explicit expo_array_commons(const std::string &s)
      {
        typedef typename Derived::value_type value_type;
        std::vector<std::string> sd;
        boost::split(sd,s,boost::is_any_of(std::string(1,derived_const_cast->separator)));
        // TODO: check here that we are not loading too many multipliers, outside expo_size_t range.
        // TODO: do it everywhere!
        const size_t w=sd.size();
        if (w == 0)
        {
          std::cout << "Warning: constructing empty expo_array." << std::endl;
          std::abort();
          return;
        }
        derived_cast->resize(w);
        for (size_t i=0;i < w;++i)
        {
          (*derived_cast)[i]=utils::lexical_converter<value_type>(sd[i]);
        }
      }
  };
}

#undef derived_const_cast
#undef derived_cast

#endif
