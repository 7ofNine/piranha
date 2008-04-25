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

#ifndef PIRANHA_POLYNOMIAL_CF_H
#define PIRANHA_POLYNOMIAL_CF_H

#include <boost/multi_index_container.hpp>

#include "../base_classes/base_series.h"
#include "../base_classes/cf_series.h"
#include "../base_classes/series_multiplication.h"
#include "../exceptions.h"
#include "../polynomial_common/common_polynomial_toolbox.h"
#include "../proxies.h"
#include "../settings.h"

#define __PIRANHA_POLYNOMIAL_CF_TP_DECL class Cf, class Expo, template <class> class I, \
  template <class, class, class, template <class> class> class Multiplier, \
  template <class> class Truncator, class Allocator
#define __PIRANHA_POLYNOMIAL_CF_TP Cf,Expo,I,Multiplier,Truncator,Allocator
#define __PIRANHA_POLYNOMIAL_CF polynomial_cf<__PIRANHA_POLYNOMIAL_CF_TP>
#define __PIRANHA_POLYNOMIAL_CF_BASE_ANCESTOR base_series<monomial<Cf,Expo,'!',Allocator>,',',Allocator,__PIRANHA_POLYNOMIAL_CF >
#define __PIRANHA_POLYNOMIAL_CF_CF_ANCESTOR cf_series< __PIRANHA_POLYNOMIAL_CF >
#define __PIRANHA_POLYNOMIAL_CF_MULT_ANCESTOR series_multiplication< __PIRANHA_POLYNOMIAL_CF, Multiplier, Truncator>
#define __PIRANHA_POLYNOMIAL_CF_COMMON_ANCESTOR common_polynomial_toolbox< __PIRANHA_POLYNOMIAL_CF >

namespace piranha
{
  template <__PIRANHA_POLYNOMIAL_CF_TP_DECL>
    class polynomial_cf:
    public __PIRANHA_POLYNOMIAL_CF_BASE_ANCESTOR,
    public __PIRANHA_POLYNOMIAL_CF_COMMON_ANCESTOR,
    public __PIRANHA_POLYNOMIAL_CF_CF_ANCESTOR,
    public __PIRANHA_POLYNOMIAL_CF_MULT_ANCESTOR
  {
      typedef monomial<Cf,Expo,'!',Allocator> term_type_;
      typedef typename term_type_::cf_type cf_type;
      typedef typename term_type_::key_type key_type;
      typedef Allocator allocator_type;
      typedef __PIRANHA_POLYNOMIAL_CF_CF_ANCESTOR cf_ancestor;
      typedef __PIRANHA_POLYNOMIAL_CF_BASE_ANCESTOR base_ancestor;
      typedef boost::multi_index_container<term_type_,typename I<term_type_>::type,allocator_type> container_type;
      typedef typename container_type::template nth_index<0>::type sorted_index;
      typedef typename container_type::template nth_index<1>::type pinpoint_index;
      friend class __PIRANHA_POLYNOMIAL_CF_CF_ANCESTOR;
      friend class __PIRANHA_POLYNOMIAL_CF_BASE_ANCESTOR;
      friend class __PIRANHA_POLYNOMIAL_CF_MULT_ANCESTOR;
      // Specify we will use the real_pow from the polynomial toolbox.
      using __PIRANHA_POLYNOMIAL_CF_COMMON_ANCESTOR::real_pow;
    public:
      // Needed typedefs.
      typedef term_type_ term_type;
      typedef typename sorted_index::const_iterator const_sorted_iterator;
      typedef typename sorted_index::iterator sorted_iterator;
      typedef typename pinpoint_index::const_iterator const_pinpoint_iterator;
      typedef typename pinpoint_index::iterator pinpoint_iterator;
      /// Default ctor.
      polynomial_cf() {nth_index<1>().max_load_factor(settings::load_factor());}
      /// Ctor from string.
      template <class ArgsTuple>
        explicit polynomial_cf(const std::string &s, const ArgsTuple &args_tuple)
      {
        nth_index<1>().max_load_factor(settings::load_factor());
        cf_ancestor::construct_from_string(s,args_tuple);
      }
      template <class ArgsTuple>
        explicit polynomial_cf(const int &n, const ArgsTuple &a)
      {
        nth_index<1>().max_load_factor(settings::load_factor());
        base_ancestor::construct_from_number(n,a);
      }
      template <class ArgsTuple>
        explicit polynomial_cf(const double &x, const ArgsTuple &a)
      {
        nth_index<1>().max_load_factor(settings::load_factor());
        base_ancestor::construct_from_number(x,a);
      }
      template <class ArgsTuple>
        explicit polynomial_cf(const psym_p &p, const int &n, const ArgsTuple &a)
      {
        nth_index<1>().max_load_factor(settings::load_factor());
        base_ancestor::construct_from_psym_p(p,n,a);
      }
      // Needed getters and setters.
      template <int N>
        typename container_type::template nth_index<N>::type &nth_index() {return m_container.template get<N>();}
      template <int N>
        const typename container_type::template nth_index<N>::type &nth_index() const {return m_container.template get<N>();}
      // TODO: place some of these methods into common polynomial toolbox?
      /// Return a vector of integers representing the polynomial.
      /**
       * If the polynomial is not a linear combination of arguments with integer coefficients, an exception will be thrown.
       * Otherwise, the polynomial's coefficients are stored into v. Used in the calculation of circular functions of 
       * Poisson series.
       */
      void get_int_linear_combination(std::vector<int> &v) const throw (unsuitable)
      {
        const const_sorted_iterator it_f = nth_index<0>().end();
        for (const_sorted_iterator it = nth_index<0>().begin(); it != it_f; ++it)
        {
          v[it->m_key.linear_arg_position()] = it->m_cf.get_int();
        }
      }
    private:
      container_type  m_container;
  };

  // Specialisation of cf mult proxy to use reference.
  template < __PIRANHA_POLYNOMIAL_CF_TP_DECL >
    class cf_mult_proxy<polynomial_cf<__PIRANHA_POLYNOMIAL_CF_TP> >:
    public reference_proxy<polynomial_cf<__PIRANHA_POLYNOMIAL_CF_TP> >
  {
      typedef reference_proxy<polynomial_cf<__PIRANHA_POLYNOMIAL_CF_TP> > ancestor;
    public:
      cf_mult_proxy():ancestor() {}
      void operator=(const polynomial_cf<__PIRANHA_POLYNOMIAL_CF_TP> &cf) {ancestor::assignment(cf);}
  };
}

#undef __PIRANHA_POLYNOMIAL_CF_TP_DECL
#undef __PIRANHA_POLYNOMIAL_CF_TP
#undef __PIRANHA_POLYNOMIAL_CF
#undef __PIRANHA_POLYNOMIAL_CF_BASE_ANCESTOR
#undef __PIRANHA_POLYNOMIAL_CF_CF_ANCESTOR
#undef __PIRANHA_POLYNOMIAL_CF_MULT_ANCESTOR
#undef __PIRANHA_POLYNOMIAL_CF_COMMON_ANCESTOR

#endif
