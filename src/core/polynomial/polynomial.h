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

#ifndef PIRANHA_POLYNOMIAL_H
#define PIRANHA_POLYNOMIAL_H

#include <boost/operators.hpp>
#include <boost/multi_index_container.hpp>
#include <memory> // For default allocator.

#include <ext/pb_ds/assoc_container.hpp>
#include <tr1/unordered_set>

#include "../base_classes/base_series.h"
#include "../base_classes/common_args_descriptions.h"
#include "../base_classes/named_series.h"
#include "../base_classes/series_multiplication.h"
#include "../polynomial_common/monomial.h"
#include "../settings_manager.h"

#define __PIRANHA_POLYNOMIAL_TP_DECL class Cf, class Expo, template <class> class I, \
  template <class, class, class, template <class, class, class> class> class Multiplier, \
  template <class, class, class> class Truncator, class Allocator
#define __PIRANHA_POLYNOMIAL_TP Cf,Expo,I,Multiplier,Truncator,Allocator
#define __PIRANHA_POLYNOMIAL polynomial<__PIRANHA_POLYNOMIAL_TP>
#define __PIRANHA_POLYNOMIAL_BASE_ANCESTOR base_series<monomial<Cf,Expo,'|',Allocator>,'\n',Allocator,__PIRANHA_POLYNOMIAL >
#define __PIRANHA_POLYNOMIAL_NAMED_ANCESTOR named_series<boost::tuple<poly_args_descr>,__PIRANHA_POLYNOMIAL >
#define __PIRANHA_POLYNOMIAL_MULT_ANCESTOR series_multiplication< __PIRANHA_POLYNOMIAL, Multiplier, Truncator>

namespace piranha
{
  template <class Term, class Allocator>
  class polynomial_container_type
  {
      typedef std::tr1::unordered_set<Term,typename Term::hasher,std::equal_to<Term>,Allocator> container_type;
//       typedef __gnu_pbds::cc_hash_table<
//         Term,
//         __gnu_pbds::null_mapped_type,
//         typename Term::hasher,
//         std::equal_to<Term>,
//         __gnu_pbds::direct_mask_range_hashing<size_t>,
//         __gnu_pbds::hash_standard_resize_policy<
//           __gnu_pbds::hash_exponential_size_policy<size_t>,
//           __gnu_pbds::hash_load_check_resize_trigger<false,size_t>,
//           true,
//           size_t
//         >,
//         true,
//         Allocator
//       > container_type;
    public:
      typedef typename container_type::const_iterator const_pinpoint_iterator;
      typedef typename container_type::iterator pinpoint_iterator;
      typedef const_pinpoint_iterator const_sorted_iterator;
      typedef pinpoint_iterator sorted_iterator;
      polynomial_container_type() {m_container.max_load_factor(.5);}
      const_pinpoint_iterator begin() const {return m_container.begin();}
      pinpoint_iterator begin() {return m_container.begin();}
      const_pinpoint_iterator end() const {return m_container.end();}
      pinpoint_iterator end() {return m_container.end();}
      void swap(polynomial_container_type &p) {m_container.swap(p.m_container);}
      size_t size() const {return m_container.size();}
      bool empty() const {return m_container.empty();}
      // TODO: check about constness below.
      pinpoint_iterator find(const Term &t) {return m_container.find(t);}
      sorted_iterator insert(const const_sorted_iterator &, const Term &t)
      {
        std::pair<pinpoint_iterator,bool> res(m_container.insert(t));
        p_assert(res.second);
        return res.first;
      }
      template <class Modifier>
        bool modify(pinpoint_iterator &it, Modifier &m)
      {
        m(*it);
        return true;
      }
      template <class Modifier>
        bool modify(const_pinpoint_iterator &it, Modifier &m)
      {
        m(*it);
        return true;
      }
      void erase(const pinpoint_iterator &it) {m_container.erase(it);}
//action_assert(derived_cast->template nth_index<0>().modify(it_new,modifier_invert_term_sign<ArgsTuple>(args_tuple)));
    private:
      container_type  m_container;
  };

  template <__PIRANHA_POLYNOMIAL_TP_DECL = std::allocator<char> >
    class polynomial:
    public __PIRANHA_POLYNOMIAL_BASE_ANCESTOR,
    public __PIRANHA_POLYNOMIAL_NAMED_ANCESTOR,
    public __PIRANHA_POLYNOMIAL_MULT_ANCESTOR,
    boost::ring_operators<__PIRANHA_POLYNOMIAL,
    boost::ring_operators<__PIRANHA_POLYNOMIAL,int,
    boost::ring_operators<__PIRANHA_POLYNOMIAL,double
    > > >
  {
      typedef monomial<Cf,Expo,'|',Allocator> term_type_;
      typedef Allocator allocator_type;
      typedef __PIRANHA_POLYNOMIAL_NAMED_ANCESTOR named_ancestor;
      typedef __PIRANHA_POLYNOMIAL_BASE_ANCESTOR base_ancestor;
//       typedef typename boost::multi_index_container <term_type_,typename I<term_type_>::type,allocator_type> container_type;
//       typedef typename container_type::template nth_index<0>::type sorted_index;
//       typedef typename container_type::template nth_index<0>::type pinpoint_index;
      typedef polynomial_container_type<term_type_,allocator_type> container_type;
      friend class __PIRANHA_POLYNOMIAL_NAMED_ANCESTOR;
      friend class __PIRANHA_POLYNOMIAL_BASE_ANCESTOR;
      friend class __PIRANHA_POLYNOMIAL_MULT_ANCESTOR;
    public:
      // Needed typedefs.
      typedef term_type_ term_type;
//       typedef typename sorted_index::const_iterator const_sorted_iterator;
//       typedef typename sorted_index::iterator sorted_iterator;
//       typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
//       typedef typename pinpoint_index::iterator pinpoint_iterator;
      typedef typename container_type::const_pinpoint_iterator const_pinpoint_iterator;
      typedef typename container_type::pinpoint_iterator pinpoint_iterator;
      typedef typename container_type::const_sorted_iterator const_sorted_iterator;
      typedef typename container_type::sorted_iterator sorted_iterator;
      // Ctors.
      polynomial() {/*nth_index<1>().max_load_factor(settings_manager::get_load_factor());*/}
      explicit polynomial(const std::string &filename)
      {
//         nth_index<1>().max_load_factor(settings_manager::get_load_factor());
        named_ancestor::construct_from_file(filename);
      }
      explicit polynomial(const int &n)
      {
//         nth_index<1>().max_load_factor(settings_manager::get_load_factor());
        base_ancestor::construct_from_number(n,named_ancestor::m_arguments);
      }
      explicit polynomial(const double &x)
      {
//         nth_index<1>().max_load_factor(settings_manager::get_load_factor());
        base_ancestor::construct_from_number(x,named_ancestor::m_arguments);
      }
      // Needed getters and setters.
      template <int>
        container_type &nth_index() {return m_container;}
      template <int>
        const container_type &nth_index() const {return m_container;}
    private:
      container_type  m_container;
  };
}

#undef __PIRANHA_POLYNOMIAL_TP_DECL
#undef __PIRANHA_POLYNOMIAL_TP
#undef __PIRANHA_POLYNOMIAL
#undef __PIRANHA_POLYNOMIAL_BASE_ANCESTOR
#undef __PIRANHA_POLYNOMIAL_NAMED_ANCESTOR
#undef __PIRANHA_POLYNOMIAL_MULT_ANCESTOR

#endif
