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

#ifndef PIRANHA_PLAIN_SERIES_MULTIPLIER_H
#define PIRANHA_PLAIN_SERIES_MULTIPLIER_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/type_traits/is_same.hpp> // For key type detection.
#include <vector>

#include "../p_assert.h"
#include "../proxies.h"
#include "../settings.h"
#include "plain_series_multiplier_mp.h"

namespace piranha
{
  /// Generic series multiplier.
  /**
   * Called from piranha::series_multiplication to multiply two piranha::base_series. This multiplier is generic,
   * provided that certain methods are implemented in the term class used for the series.
   *
   * This class can be extended to build more specific multipliers.
   */
  template <class Series1, class Series2, class ArgsTuple, template <class> class Truncator>
    class plain_series_multiplier
  {
      friend struct base_insert_multiplication_result;
      friend class Truncator<plain_series_multiplier>;
    public:
      // These typedefs are public because truncators may want to use them.
      /// Alias for term type of first input series and return value series.
      typedef typename Series1::term_type term_type1;
      /// Alias for term type of second input series.
      typedef typename Series2::term_type term_type2;
      /// Alias for the truncator type.
      typedef Truncator<plain_series_multiplier> truncator_type;
    private:
      typedef boost::multi_index_container
      <
        term_type1,
        boost::multi_index::indexed_by
        <
          boost::multi_index::hashed_unique<boost::multi_index::identity<term_type1> >
        >
      >
      mult_set;
      BOOST_STATIC_ASSERT((boost::is_same<typename term_type1::key_type,typename term_type2::key_type>::value));
      typedef cf_mult_proxy<typename term_type1::cf_type> cf_proxy_type1;
      typedef cf_mult_proxy<typename term_type2::cf_type> cf_proxy_type2;
      typedef key_mult_proxy<typename term_type1::key_type> key_proxy_type;
    public:
      plain_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        m_s1(s1),m_s2(s2),m_args_tuple(args_tuple),m_size1(m_s1.template nth_index<0>().size()),
        m_size2(m_s2.template nth_index<0>().size()),m_retval(retval),m_cfs1(m_size1),m_cfs2(m_size2),
        m_keys1(m_size1),m_keys2(m_size2),m_trunc(*this)
      {
        // Set proper load factor for hash set.
        m_set.max_load_factor(settings::load_factor());
        // Cache the series terms into separate vectors for coefficients and keys.
        cache_series_terms(m_s1,m_cfs1,m_keys1);
        cache_series_terms(m_s2,m_cfs2,m_keys2);
      }
      /// Perform multiplication.
      /**
       * Method called by piranha::series_multiplication. Can be overloaded by derived multipliers.
       * Internally it just calls perform_plain_multiplication.
       */
      void perform_multiplication()
      {
        perform_plain_multiplication();
      }
    protected:
      void perform_plain_multiplication()
      {
        plain_multiplication();
        plain_insert_result_into_retval();
      }
    private:
      // Perform plain multiplication.
      void plain_multiplication()
      {
        typedef typename term_type1::multiplication_result mult_res;
        mult_res res;
        for (size_t i = 0; i < m_size1; ++i)
        {
          for (size_t j = 0; j < m_size2; ++j)
          {
            if (m_trunc.skip(m_cfs1[i].get(),m_keys1[i].get(),m_cfs2[j].get(),m_keys2[j].get(),*this))
            {
              break;
            }
            term_type1::multiply(m_cfs1[i].get(),m_keys1[i].get(),m_cfs2[j].get(),m_keys2[j].get(),res,m_args_tuple);
            insert_multiplication_result<mult_res>::run(res,*this);
          }
        }
      }
      template <class Series>
        void cache_series_terms(const Series &s,
        std::vector<cf_mult_proxy<typename Series::term_type::cf_type> > &cfs,
        std::vector<key_mult_proxy<typename Series::term_type::key_type> > &keys)
      {
        typedef typename Series::const_sorted_iterator const_sorted_iterator;
        const const_sorted_iterator it_f = s.template nth_index<0>().end();
        size_t i=0;
        for (const_sorted_iterator it = s.template nth_index<0>().begin(); it != it_f; ++it)
        {
          cfs[i] = it->m_cf;
          keys[i] = it->m_key;
          ++i;
        }
      }
      // After the multiplication has been performed and the result stored in the temporary hash table,
      // fetch the terms from there and put them into retval.
      void plain_insert_result_into_retval()
      {
        typedef typename mult_set::const_iterator hash_iterator;
        typedef typename Series1::sorted_iterator sorted_iterator;
        term_type1 term;
        sorted_iterator it_hint = m_retval.template nth_index<0>().end();
        const hash_iterator it_f = m_set.end();
        for (hash_iterator it = m_set.begin(); it != it_f; ++it)
        {
          term.m_cf = it->m_cf;
          term.m_key = it->m_key;
          it_hint = m_retval.template insert<false,true>(term,m_args_tuple,it_hint);
        }
      }
    protected:
      // References to the series.
      const Series1                 &m_s1;
      const Series2                 &m_s2;
      // Reference to the arguments tuple.
      const ArgsTuple               &m_args_tuple;
      // Sizes of the series.
      const size_t                  m_size1;
      const size_t                  m_size2;
      Series1                       &m_retval;
      // Vectors of input coefficients converted for representation during series multiplication.
      std::vector<cf_proxy_type1>   m_cfs1;
      std::vector<cf_proxy_type2>   m_cfs2;
      // Vectors of input keys converted for representation during series multiplication.
      std::vector<key_proxy_type>   m_keys1;
      std::vector<key_proxy_type>   m_keys2;
      // Container to store the result of the multiplications.
      mult_set                      m_set;
      // Truncator. This must be the last one defined because it will take *this
      // as parameter for construction.
      truncator_type                m_trunc;
  };
}

#endif
