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

#ifndef PIRANHA_POISSON_SERIES_H
#define PIRANHA_POISSON_SERIES_H

#include <boost/multi_index_container.hpp>
#include <memory> // For default allocator.

#include "../arg_manager.h"
#include "../base_classes/base_series.h"
#include "../base_classes/named_series.h"
#include "../poisson_series_common/poisson_series_term.h"
#include "../polynomial_cf/polynomial_cf.h"

#define __PIRANHA_POISSON_SERIES_TP_DECL class Cf, class Expo, class Trig, template <class> class IPoly,\
  template <class> class ITrig, class Allocator
#define __PIRANHA_POISSON_SERIES_TP Cf,Expo,Trig,IPoly,ITrig,Allocator
#define __PIRANHA_POISSON_SERIES poisson_series<__PIRANHA_POISSON_SERIES_TP>
#define __PIRANHA_POISSON_SERIES_BASE_ANCESTOR base_series<poisson_series_term<polynomial_cf<Cf,Expo,IPoly,Allocator>,Trig,'|',Allocator>,\
  '\n',Allocator,__PIRANHA_POISSON_SERIES > 
#define __PIRANHA_POISSON_SERIES_NAMED_ANCESTOR named_series<boost::tuple<poly_args_descr,trig_args_descr>,__PIRANHA_POISSON_SERIES >

namespace piranha
{
  template <__PIRANHA_POISSON_SERIES_TP_DECL = std::allocator<char> >
    class poisson_series:
    protected __PIRANHA_POISSON_SERIES_BASE_ANCESTOR,
    public __PIRANHA_POISSON_SERIES_NAMED_ANCESTOR
  {
      typedef poisson_series_term<polynomial_cf<Cf,Expo,IPoly,Allocator>,Trig,'|',Allocator> term_type_;
      typedef Allocator allocator_type;
      typedef __PIRANHA_POISSON_SERIES_NAMED_ANCESTOR named_ancestor;
      typedef __PIRANHA_POISSON_SERIES_BASE_ANCESTOR base_ancestor;
      typedef typename boost::multi_index_container <term_type_,
        typename ITrig<term_type_>::type,allocator_type> container_type;
      typedef typename container_type::template nth_index<0>::type sorted_index;
      typedef typename container_type::template nth_index<1>::type pinpoint_index;
      friend class __PIRANHA_POISSON_SERIES_NAMED_ANCESTOR;
      friend class __PIRANHA_POISSON_SERIES_BASE_ANCESTOR;
    public:
      // Needed typedefs.
      typedef term_type_ term_type;
      typedef typename sorted_index::const_iterator const_sorted_iterator;
      typedef typename sorted_index::iterator sorted_iterator;
      typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
      typedef typename pinpoint_index::iterator pinpoint_iterator;
      // Ctors.
      poisson_series() {}
      explicit poisson_series(const std::string &filename) {named_ancestor::construct_from_file(filename);}
      explicit poisson_series(const int &n) {base_ancestor::construct_from_number(n,named_ancestor::m_arguments);}
      explicit poisson_series(const double &x) {base_ancestor::construct_from_number(x,named_ancestor::m_arguments);}
      // Needed getters and setters.
      template <int N>
        typename container_type::template nth_index<N>::type &nth_index() {return m_container.template get<N>();}
      template <int N>
        const typename container_type::template nth_index<N>::type &nth_index() const {return m_container.template get<N>();}
    private:
      // Private ctors.
      explicit poisson_series(const typename named_ancestor::args_tuple_type &args_tuple) {named_ancestor::construct_from_args(args_tuple);}
    private:
      container_type  m_container;
  };
}

#undef __PIRANHA_POISSON_SERIES_TP_DECL
#undef __PIRANHA_POISSON_SERIES_TP
#undef __PIRANHA_POISSON_SERIES
#undef __PIRANHA_POISSON_SERIES_BASE_ANCESTOR
#undef __PIRANHA_POISSON_SERIES_NAMED_ANCESTOR

#endif
