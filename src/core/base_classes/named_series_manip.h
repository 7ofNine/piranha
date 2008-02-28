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

  // Meta-programming for appending an argument.
  template <class ArgsDescr>
    struct named_series_append_arg
  {
    static void run(const std::string &s,
      typename ntuple<vector_psym_p,boost::tuples::length<ArgsDescr>::value>::type &args_tuple,
      const psym_p &arg)
    {
      switch (ArgsDescr::head_type::name == s)
      {
        case true:
          // Check that the argument is not already present in this set.
          for (vector_psym_p::iterator it = args_tuple.get_head().begin(); it != args_tuple.get_head().end(); ++it)
          {
            if (arg == (*it))
            {
              std::cout << "Error: " << s << " argument '" << (*it)->name() << "' already present in the set." << std::endl;
              return;
            }
          }
          args_tuple.get_head().push_back(arg);
          break;
        case false:
          named_series_append_arg<typename ArgsDescr::tail_type>::run(s,args_tuple.get_tail(),arg);
      }
    }
  };

  template <>
    struct named_series_append_arg<boost::tuples::null_type>
  {
    static void run(const std::string &s, const boost::tuples::null_type &, const psym_p &)
    {
      std::cout << "Error: '" << s << "' arguments are not known." << std::endl;
    }
  };

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::append_arg(const std::string &s, const psym_p &arg)
  {
    p_assert(derived_const_cast->empty());
    named_series_append_arg<arguments_description>::run(s,m_arguments,arg);
  }
}

#endif
