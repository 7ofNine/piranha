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

#ifndef PIRANHA_NAMED_SERIES_MANIP_H
#define PIRANHA_NAMED_SERIES_MANIP_H

#include <boost/static_assert.hpp>

namespace piranha
{
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <bool AdditionalChecks, bool Sign, class Term2, class SortedIterator>
    inline SortedIterator named_series<__PIRANHA_NAMED_SERIES_TP>::insert(const Term2 &term, SortedIterator it_hint)
  {
    return derived_cast->insert<AdditionalChecks,Sign,Term2,arguments_tuple_type>(term,m_arguments,it_hint);
  }

  /// Perform plain insertion with all checks.
  /**
   * Simple wrapper around named_series::insert.
   */
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Term2, class SortedIterator>
    inline SortedIterator
    named_series<__PIRANHA_NAMED_SERIES_TP>::insert(const Term2 &term, SortedIterator it_hint)
  {
    return insert<true,true>(term,m_arguments,it_hint);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <int N, class Iterator>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::term_erase(Iterator it)
  {
    derived_cast->template term_erase<N>(m_arguments,it);
  }

  // Meta-programming for appending argument.
  template <int N>
    struct append_series_argument
  {
    BOOST_STATIC_ASSERT(N > 0);
    template <class ArgsTuple, class ArgsDescr>
      static void run(const std::string &s, ArgsTuple &args_tuple, const psym_p &arg)
    {
      switch (boost::tuples::element<N,ArgsDescr>::type::name == s)
      {
        case true:
          args_tuple.template get<N>().push_back(arg);
          break;
        case false:
          append_series_argument<N-1>::template run<ArgsTuple,ArgsDescr>(s,args_tuple,arg);
      }
    }
  };

  template <>
    struct append_series_argument<0>
  {
    template <class ArgsTuple, class ArgsDescr>
      static void run(const std::string &s, ArgsTuple &args_tuple, const psym_p &arg)
    {
      switch (boost::tuples::element<0,ArgsDescr>::type::name == s)
      {
        case true:
          args_tuple.template get<0>().push_back(arg);
          break;
        case false:
          std::cout << "Error: '" << s << "' arguments are not known." << std::endl;
      }
    }
  };

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::append_arg(const std::string &s, const psym_p &arg)
  {
    p_assert(derived_const_cast->empty());
    append_series_argument<n_arguments_sets-1>::
      template run<arguments_tuple_type,arguments_description>(s,m_arguments,arg);
  }
}

#endif
