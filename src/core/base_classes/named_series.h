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

#ifndef PIRANHA_NAMED_SERIES_H
#define PIRANHA_NAMED_SERIES_H

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>

#include "../ntuple.h"
#include "../psymbol.h"

// Useful shortcuts.
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_NAMED_SERIES_TP_DECL class ArgsDescr, class Derived
#define __PIRANHA_NAMED_SERIES_TP ArgsDescr,Derived

namespace piranha
{
  /// Named series toolbox.
  /**
   * Toolbox for generating series with arguments.
   * ArgsDescr must be a boost::tuple of structures each one containing a static const string
   * called "name" naming the arguments of the series.
   */
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    class named_series
  {
      typedef ArgsDescr arguments_description;
    public:
      /// Compile-time constant for the number of arguments sets.
      static const int n_arguments_sets = boost::tuples::length<arguments_description>::value;
      typedef ntuple<vector_psym_p,n_arguments_sets> arguments_tuple_type;
      BOOST_STATIC_ASSERT(n_arguments_sets > 0);
      template <bool, bool, class Term2, class SortedIterator>
        SortedIterator insert(const Term2 &, SortedIterator);
      template <class Term2, class SortedIterator>
        SortedIterator insert(const Term2 &, SortedIterator);
      template <int N, class Iterator>
        void term_erase(Iterator);
    //private:
      static int args_pos(const std::string &);
    protected:
      arguments_tuple_type  m_arguments;
  };

  // Initialization of static members.
   template <__PIRANHA_NAMED_SERIES_TP_DECL>
     const int named_series<__PIRANHA_NAMED_SERIES_TP>::n_arguments_sets;
}

#include "named_series_manip.h"
#include "named_series_probe.h"

#undef derived_const_cast
#undef derived_cast
#undef __PIRANHA_NAMED_SERIES_TP_DECL
#undef __PIRANHA_NAMED_SERIES_TP

#endif
