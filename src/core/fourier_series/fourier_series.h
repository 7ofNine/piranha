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

#ifndef PIRANHA_FOURIER_SERIES_H
#define PIRANHA_FOURIER_SERIES_H

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>
#include <memory> // For default allocator.

#include "../arg_manager.h"
#include "../base_classes/base_series.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/named_series.h"
#include "../fourier_series_common/fs_term.h"
#include "../ntuple.h"

namespace piranha
{
  /// Norm extractor.
  /**
   * Extract norm from term, asserting that arguments have been assigned through piranha::arg_manager.
   */
  template <class Term, class ArgsTuple>
    struct norm_extractor
  {
    typedef double result_type;
    double operator()(const Term &t) const
    {
      p_assert(arg_manager<ArgsTuple>::assigned());
      return t.m_cf.norm(arg_manager<ArgsTuple>::get());
    }
  };

  /// Norm-based indices for base_pseries.
  /**
   * This class specifies the following indices to be used in piranha::base_pseries: a hashed index for the
   * identification
   * of terms and a norm-sorted index to discard terms in multiplications. The class is to be used as the I
   * parameter in piranha::base_pseries classes.
   */
  template <class Term>
    struct fourier_series_norm_index
  {
    typedef typename Term::key_type trig_type;
    typedef boost::multi_index::indexed_by <
      boost::multi_index::ordered_unique <
      boost::multi_index::composite_key <
      Term,
      norm_extractor<Term,typename ntuple<vector_psym_p,1>::type>,
      key_extractor<Term>
      >,
      boost::multi_index::composite_key_compare<
      std::greater<double>,
      std::less<trig_type>
      >
      >,
      boost::multi_index::hashed_unique<boost::multi_index::identity<Term> >
      > type;
  };


#define __PIRANHA_FOURIER_SERIES_TP_DECL class Cf, class Trig, template <class> class I, class Allocator
#define __PIRANHA_FOURIER_SERIES_TP Cf,Trig,I,Allocator
#define __PIRANHA_FOURIER_SERIES fourier_series<__PIRANHA_FOURIER_SERIES_TP>
#define __PIRANHA_FOURIER_SERIES_BASE_ANCESTOR base_series<fs_term<Cf,Trig,'|',Allocator>,'\n',Allocator,__PIRANHA_FOURIER_SERIES >
#define __PIRANHA_FOURIER_SERIES_NAMED_ANCESTOR named_series<boost::tuple<trig_args_descr>,__PIRANHA_FOURIER_SERIES >

  template <__PIRANHA_FOURIER_SERIES_TP_DECL = std::allocator<char> >
    class fourier_series:
    protected __PIRANHA_FOURIER_SERIES_BASE_ANCESTOR,
    public __PIRANHA_FOURIER_SERIES_NAMED_ANCESTOR
  {
      typedef fs_term<Cf,Trig,'|',Allocator> term_type_;
      typedef Allocator allocator_type;
      typedef __PIRANHA_FOURIER_SERIES_NAMED_ANCESTOR named_ancestor;
      typedef __PIRANHA_FOURIER_SERIES_BASE_ANCESTOR base_ancestor;
      typedef typename boost::multi_index_container <term_type_,typename I<term_type_>::type,allocator_type> container_type;
      typedef typename container_type::template nth_index<0>::type sorted_index;
      typedef typename container_type::template nth_index<1>::type pinpoint_index;
      friend class __PIRANHA_FOURIER_SERIES_NAMED_ANCESTOR;
      friend class __PIRANHA_FOURIER_SERIES_BASE_ANCESTOR;
    public:
      // Needed typedefs.
      typedef term_type_ term_type;
      typedef typename sorted_index::const_iterator const_sorted_iterator;
      typedef typename sorted_index::iterator sorted_iterator;
      typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
      typedef typename pinpoint_index::iterator pinpoint_iterator;
      // Ctors.
      fourier_series() {}
      explicit fourier_series(const std::string &filename)
      {
        named_ancestor::read_from_file(filename);
      }
      // Needed getters and setters.
      template <int N>
        typename container_type::template nth_index<N>::type &nth_index() {return m_container.template get<N>();}
      template <int N>
        const typename container_type::template nth_index<N>::type &nth_index() const {return m_container.template get<N>();}
    private:
      container_type  m_container;
  };

#undef __PIRANHA_FOURIER_SERIES_TP_DECL
#undef __PIRANHA_FOURIER_SERIES_TP
#undef __PIRANHA_FOURIER_SERIES
#undef __PIRANHA_FOURIER_SERIES_BASE_ANCESTOR
#undef __PIRANHA_FOURIER_SERIES_NAMED_ANCESTOR
}

#endif
