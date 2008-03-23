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

#ifndef PIRANHA_SERIES_MULTIPLIER_H
#define PIRANHA_SERIES_MULTIPLIER_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <valarray>

#include "../settings_manager.h"
#include "../type_traits.h"
#include "series_multiplier_mp.h"

namespace piranha
{
  template <class Series1, class Series2, class ArgsTuple, template <class, class, class> class Truncator>
    class series_multiplier
  {
      typedef Truncator<Series1,Series2,ArgsTuple> truncator_type;
    protected:
      typedef typename Series1::term_type term_type;
      typedef typename Series1::term_type::cf_type cf1;
      typedef typename Series2::term_type::cf_type cf2;
      typedef typename Series1::term_type::key_type key;
      typedef typename constant_copy<cf1>::type sm_cf1;
      typedef typename constant_copy<cf2>::type sm_cf2;
      typedef typename constant_copy<key>::type sm_key;
      typedef boost::multi_index_container
      <
        term_type,
        boost::multi_index::indexed_by
        <
          boost::multi_index::hashed_unique<boost::multi_index::identity<term_type> >
        >
      >
      mult_set;
    public:
      series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        m_s1(s1),m_s2(s2),m_retval(retval),m_args_tuple(args_tuple),m_cfs1(),m_cfs2(),m_keys1(),m_keys2(),m_set()
      {
        // Set proper load factor for hash set.
        m_set.max_load_factor(settings_manager::get_load_factor());
      }
      void perform_multiplication()
      {
        adjust_sizes();
        cache_series_terms(m_s1,m_cfs1,m_keys1);
        cache_series_terms(m_s2,m_cfs2,m_keys2);
        plain_multiplication(m_args_tuple);
        plain_insert_result_into_retval(m_args_tuple);
      }
    protected:
      template <class Series>
        void cache_series_terms(const Series &s,
        std::valarray<typename constant_copy<typename Series::term_type::cf_type>::type> &cfs,
        std::valarray<typename constant_copy<typename Series::term_type::key_type>::type> &keys)
      {
        typedef typename Series::const_sorted_iterator const_sorted_iterator;
        const const_sorted_iterator it_f = s.template nth_index<0>().end();
        size_t i=0;
        for (const_sorted_iterator it = s.template nth_index<0>().begin(); it != it_f; ++it)
        {
          constant_copy<typename Series::term_type::cf_type>::assign(cfs[i],it->m_cf); 
          constant_copy<typename Series::term_type::key_type>::assign(keys[i],it->m_key);
          ++i;
        }
      }
      // Perform plain multiplication.
      void plain_multiplication(const ArgsTuple &args_tuple)
      {
        typedef typename term_type::multiplication_result mult_res;
        truncator_type trunc(m_s1,m_s2,args_tuple);
        const size_t size1 = m_cfs1.size(), size2 = m_cfs2.size();
        p_assert(size1 == m_keys1.size() and size2 == m_keys2.size());
        mult_res res;
        for (size_t i = 0; i < size1; ++i)
        {
          for (size_t j = 0; j < size2; ++j)
          {
            if (trunc.skip_from_here(
              constant_copy<cf1>::get(m_cfs1[i]),
              constant_copy<key>::get(m_keys1[i]),
              constant_copy<cf2>::get(m_cfs2[j]),
              constant_copy<key>::get(m_keys2[j])))
            {
              break;
            }
            term_type::multiply(
              constant_copy<cf1>::get(m_cfs1[i]),
              constant_copy<key>::get(m_keys1[i]),
              constant_copy<cf2>::get(m_cfs2[j]),
              constant_copy<key>::get(m_keys2[j]),
              res,args_tuple);
            insert_multiplication_result<mult_res>::run(res,m_set,args_tuple,trunc);
          }
        }
      }
      // After the multiplication has been performed and the result stored in the temporary hash table,
      // fetch the terms from there and put them into retval.
      void plain_insert_result_into_retval(const ArgsTuple &args_tuple)
      {
        typedef typename mult_set::const_iterator hash_iterator;
        typedef typename Series1::const_sorted_iterator sorted_iterator;
        term_type term;
        sorted_iterator it_hint = m_retval.template nth_index<0>().end();
        const hash_iterator it_f = m_set.end();
        for (hash_iterator it = m_set.begin(); it != it_f; ++it)
        {
          term.m_cf = it->m_cf;
          term.m_key = it->m_key;
          it_hint = m_retval.template insert<false,true>(term,args_tuple,it_hint);
        }
      }
      void adjust_sizes()
      {
        const size_t size1 = m_s1.template nth_index<0>().size(), size2 = m_s2.template nth_index<0>().size();
        m_cfs1.resize(size1);
        m_cfs2.resize(size2);
        m_keys1.resize(size1);
        m_keys2.resize(size2);
      }
    protected:
      // References to the series.
      const Series1           &m_s1;
      const Series2           &m_s2;
      Series1                 &m_retval;
      const ArgsTuple         &m_args_tuple;
      // Vectors of input coefficients converted for representation during series multiplication.
      std::valarray<sm_cf1>   m_cfs1;
      std::valarray<sm_cf2>   m_cfs2;
      // Vectors of input keys converted for representation during series multiplication.
      std::valarray<sm_key>   m_keys1;
      std::valarray<sm_key>   m_keys2;
      // Container to store the result of the multiplications.
      mult_set                m_set;
  };
}

#endif
