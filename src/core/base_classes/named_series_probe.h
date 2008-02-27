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

#ifndef PIRANHA_NAMED_SERIES_PROBE_H
#define PIRANHA_NAMED_SERIES_PROBE_H

namespace piranha
{
  // Meta-programmed structure for the identification of arguments' names.
  template <int N>
    struct arguments_identifier
  {
    BOOST_STATIC_ASSERT(N > 0);
    template <class ArgsDescr>
      static int run(const std::string &s)
    {
      switch (boost::tuples::element<N,ArgsDescr>::type::name == s)
      {
        case true:
          return N;
        case false:
          return arguments_identifier<N-1>::template run<ArgsDescr>(s);
      };
    }
  };

  template <>
    struct arguments_identifier<0>
  {
    template <class ArgsDescr>
      static int run(const std::string &s)
    {
      switch (boost::tuples::element<0,ArgsDescr>::type::name == s)
      {
        case true:
          return 0;
        default:
          return -1;
      }
    }
  };

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline int named_series<__PIRANHA_NAMED_SERIES_TP>::args_pos(const std::string &s)
  {
    return arguments_identifier<n_arguments_sets-1>::template run<arguments_description>(s);
  }
}

#endif
