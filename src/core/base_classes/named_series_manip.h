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

#include "../config.h" // For (un)likely

namespace piranha
{
  // Template metaprogramming for swapping of arguments in named series.
  template <class ArgsTuple>
    struct named_series_swap
  {
    static void run(ArgsTuple &a1, ArgsTuple &a2)
    {
      a1.get_head().swap(a2.get_head());
      named_series_swap<typename ArgsTuple::tail_type>::run(a1.get_tail(),
        a2.get_tail());
    }
  };

  template <>
    struct named_series_swap<boost::tuples::null_type>
  {
    static void run(const boost::tuples::null_type &, const boost::tuples::null_type &)
    {}
  };

  /// Swap contents of series.
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::swap(Derived &ps2)
  {
    // Do something only if we are not swapping with self.
    if (derived_cast != &ps2)
    {
      named_series_swap<args_tuple_type>::run(m_arguments,ps2.m_arguments);
      derived_cast->swap_terms(ps2);
    }
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
    p_assert(derived_const_cast->template nth_index<0>().empty());
    named_series_append_arg<arguments_description>::run(s,m_arguments,arg);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::merge_args(const Derived2 &ps2)
  {
    if (static_cast<void *>(this) == static_cast<void const *>(&ps2))
    {
      std::cout << "Trying to merge with self, returning." << std::endl;
      return;
    }
    if (unlikely(!is_args_compatible(ps2)))
    {
      merge_incompatible_args(ps2);
    }
  }

  // Template metaprogramming for getting the relative layout of two series.
  template <class ArgsTuple>
    struct named_series_get_layout
  {
    static void run(const ArgsTuple &a1, const ArgsTuple &a2,
      typename ntuple<std::vector<std::pair<bool,size_t> >,
      boost::tuples::length<ArgsTuple>::value>::type &l)
    {
      const size_t size1 = a1.get_head().size(), size2 = a2.get_head().size();
      // First we must construct a2's layout wrt to a1.
      l.get_head().resize(size2);
      for (size_t i=0;i < size2;++i)
      {
        // If we won't find a2's element, we'll mark it as not found.
        l.get_head()[i].first=false;
        // For each of a2's elements, look for that same element in a1.
        for (size_t j=0;j < size1;++j)
        {
          if (a1.get_head()[j] == a2.get_head()[i])
          {
            // We found it, mark as found and proceed to next a2 element.
            l.get_head()[i].first=true;
            l.get_head()[i].second=j;
            break;
          }
        }
      }
      // Now we must take care of those elements of a1 that are not represented in the layout (i.e., they are not in a2)
      for (size_t i=0;i < size1;++i)
      {
        // Look for element index i in the layout.
        bool found = false;
        const size_t l_size = l.get_head().size();
        for (size_t j=0;j < l_size;++j)
        {
          if (l.get_head()[j].first and l.get_head()[j].second == i)
          {
            found = true;
            break;
          }
        }
        // If we did not find it, append it to the layout.
        if (!found)
        {
          l.get_head().push_back(std::pair<bool,size_t>(true,i));
        }
      }
      named_series_get_layout<typename ArgsTuple::tail_type>::run(a1.get_tail(),
        a2.get_tail(),l.get_tail());
    }
  };

  template <>
    struct named_series_get_layout<boost::tuples::null_type>
  {
    static void run(const boost::tuples::null_type &, const boost::tuples::null_type &,
      const boost::tuples::null_type &)
    {}
  };

  // Template metaprogramming for applying a layout to a series.
  template <class ArgsTuple>
    struct named_series_apply_layout_to_args
  {
    static void run(ArgsTuple &a1, const ArgsTuple &a2,
      const typename ntuple<std::vector<std::pair<bool,size_t> >,
      boost::tuples::length<ArgsTuple>::value>::type &l)
    {
      const size_t l_size = l.get_head().size();
      // The layout must have at least all arguments in v1.
      p_assert(l_size >= a1.get_head().size());
      // Memorize the old vector.
      const vector_psym_p old(a1.get_head());
      // Make space.
      a1.get_head().resize(l_size);
      for (size_t i=0; i < l_size; ++i)
      {
        switch (l.get_head()[i].first)
        {
          case true:
            // The argument was present in the old arguments sets. Copy it over.
            p_assert(l.get_head()[i].second < old.size());
            a1.get_head()[i] = old[l.get_head()[i].second];
            break;
          case false:
            // The argument was not present in the old arguments sets. Fetch it from a2.
            p_assert(i < a2.get_head().size());
            a1.get_head()[i] = a2.get_head()[i];
        }
      }
      named_series_apply_layout_to_args<typename ArgsTuple::tail_type>::run(a1.get_tail(),
        a2.get_tail(),l.get_tail());
    }
  };

  template <>
    struct named_series_apply_layout_to_args<boost::tuples::null_type>
  {
    static void run(const boost::tuples::null_type &, const boost::tuples::null_type &,
      const boost::tuples::null_type &)
    {}
  };

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    template <class Derived2>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::merge_incompatible_args(const Derived2 &ps2)
  {
    // Build an empty retval and assign it the same arguments as this.
    Derived retval;
    retval.m_arguments = m_arguments;
    // Build a tuple of layouts.
    typename ntuple<std::vector<std::pair<bool,size_t> >,n_arguments_sets>::type l;
    // Get the relative layouts of this wrt ps2 and put the result into l.
    named_series_get_layout<args_tuple_type>::run(retval.m_arguments,ps2.m_arguments,l);
    // Apply the layout to the arguments tuple of retval.
    named_series_apply_layout_to_args<args_tuple_type>::run(retval.m_arguments,ps2.m_arguments,l);
    // Apply the layout to all terms of this, which will be inserted into retval.
    derived_cast->apply_layout_to_terms(retval.m_arguments,l,retval);
    // Finally, swap the contents of retval with this.
    swap(retval);
  }
}

#endif
