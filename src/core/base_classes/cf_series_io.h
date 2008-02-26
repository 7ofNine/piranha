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

#ifndef PIRANHA_CF_SERIES_IO_H
#define PIRANHA_CF_SERIES_IO_H

namespace piranha
{
  template <__PIRANHA_CF_SERIES_TP_DECL>
    template <class ArgsTuple>
    inline cf_series<__PIRANHA_CF_SERIES_TP>::cf_series(std::string &str, const ArgsTuple &args_tuple)
  {
    typedef typename Derived::term_type term_type;
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    const char separator = Derived::separator;
    // Remove extra spaces.
    boost::trim(str);
    // Split into single terms.
    std::vector<std::string> vs;
    boost::split(vs,str,boost::is_any_of(std::string(separator)));
    // Fetch the end of the sorted index as hint.
    const_sorted_iterator it_hint = derived_const_cast->g_s_index().end();
    const size_t length = vs.size();
    for (size_t i=0; i < length; ++i)
    {
      try
      {
        // Try to build the term from the string.
        term_type term(vs[i],args_tuple);
        it_hint = derived_cast->insert(term,args_tuple,it_hint);
      }
      catch(bad_input &b)
      {
        std::cout << b.what() << std::endl;
        std::cout << "Error reading term in coefficient series, ignoring it." << std::endl;
      }
    }
  }
}

#endif
