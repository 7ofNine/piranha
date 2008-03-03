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

#ifndef PIRANHA_BASE_SERIES_IO_H
#define PIRANHA_BASE_SERIES_IO_H

#include "../stream_manager.h"

namespace piranha
{
  template <__PIRANHA_BASE_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline void base_series<__PIRANHA_BASE_SERIES_TP>::print_terms_plain(std::ostream &stream,
    const ArgsTuple &args_tuple, int limit) const
  {
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    stream_manager::setup_print(stream);
    size_t j=0, lim;
    if (limit < 0)
    {
      lim=derived_const_cast->template nth_index<0>().size();
    }
    else
    {
      lim=(size_t)limit;
    }
    const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
    for (const_sorted_iterator it=derived_const_cast->template nth_index<0>().begin();it!=it_f;++it)
    {
      if (j == lim)
      {
        break;
      }
      it->print_plain(stream,args_tuple);
      stream << separator;
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
