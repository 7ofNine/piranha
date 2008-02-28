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

#ifndef PIRANHA_NAMED_SERIES_IO_H
#define PIRANHA_NAMED_SERIES_IO_H

#include <boost/tuple/tuple.hpp>
#include <iostream>

namespace piranha
{
  template <class ArgsDescr>
    inline void named_series_print_plain(std::ostream &stream,
    const typename ntuple<vector_psym_p,boost::tuples::length<ArgsDescr>::value>::type &args_tuple)
  {
    for (size_t i=0; i < args_tuple.get_head().size(); ++i)
    {
      stream << "[" << ArgsDescr::head_type::name << "_arg]" << std::endl;
      args_tuple.get_head()[i]->print(stream);
    }
    named_series_print_plain<typename ArgsDescr::tail_type>(stream,args_tuple.get_tail());
  }

  template <>
    inline void named_series_print_plain<boost::tuples::null_type>(std::ostream &,
    const ntuple<vector_psym_p,boost::tuples::length<boost::tuples::null_type>::value>::type &)
  {}

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print_plain(std::ostream &stream) const
  {
    named_series_print_plain<arguments_description>(stream,m_arguments);
  }
}

#endif
