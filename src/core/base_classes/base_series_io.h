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

#ifndef PIRANHA_BASE_SERIES_H
#define PIRANHA_BASE_SERIES_H

#include "../stream_manager.h"

namespace piranha
{
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms(std::ostream &stream,
    const ArgsTuple &args_tuple, int limit) const
  {
    switch (stream_manager::format())
    {
      case stream_manager::plain:
        print_terms_plain(std::cout,limit);
        break;
      case stream_manager::latex:
        print_terms_latex(std::cout,limit);
    }
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_plain(std::ostream &stream,
    const ArgsTuple &args_tuple, int limit) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    stream_manager::setup_print(out_stream);
    size_t j=0, lim;
    if (n < 0)
    {
      lim=derived_const_cast->length();
    }
    else
    {
      lim=(size_t)n;
    }
    const const_sorted_iterator it_f = derived_const_cast->g_s_index().end();
    for (const_sorted_iterator it=derived_const_cast->g_s_index().begin();it!=it_f;++it)
    {
      if (j == lim)
      {
        break;
      }
      it->print_plain(out_stream,args_tuple);
      out_stream << separator;
      ++j;
    }
  }

  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_latex(std::ostream &stream,
    const ArgsTuple &args_tuple, int limit) const
  {
// TODO: to be implemented.
  }
}

#endif
