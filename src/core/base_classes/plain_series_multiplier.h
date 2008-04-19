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
#include <valarray>

#include "../p_assert.h"
#include "../settings.h"
#include "../type_traits.h"
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
    public:
      // These typedefs are public because truncators may want to use them.
      /// Alias for term type of first input series and return value series.
      typedef typename Series1::term_type term_type;
      /// Alias for the coefficient type of the first input series.
      typedef typename Series1::term_type::cf_type cf_type1;
      /// Alias for the coefficient type of the second input series.
      typedef typename Series2::term_type::cf_type cf_type2;
      /// Alias for the key type (common to both input series).
      typedef typename Series1::term_type::key_type key_type;
      /// Alias for the truncator type.
      typedef Truncator<plain_series_multiplier> truncator_type;
      typedef boost::multi_index_container
      <
        term_type,
        boost::multi_index::indexed_by
        <
          boost::multi_index::hashed_unique<boost::multi_index::identity<term_type> >
        >
      >
      mult_set;
    private:
      typedef typename series_mult_rep<cf_type1>::type sm_cf1;
      typedef typename series_mult_rep<cf_type2>::type sm_cf2;
      typedef typename series_mult_rep<key_type>::type sm_key;
    public:
      plain_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        m_s1(s1),m_s2(s2),m_args_tuple(args_tuple),m_size1(m_s1.template nth_index<0>().size()),
        m_size2(m_s2.template nth_index<0>().size()),m_retval(retval),m_trunc(*this)
      {
        // Set proper load factor for hash set.
        m_set.max_load_factor(settings::load_factor());
      }
      void perform_multiplication()
      {
        perform_plain_multiplication();
      }
    protected:
      void perform_plain_multiplication()
      {
        adjust_input_sizes();
        cache_series_terms(m_s1,m_cfs1,m_keys1);
        cache_series_terms(m_s2,m_cfs2,m_keys2);
        plain_multiplication();
        plain_insert_result_into_retval();
      }
      template <class Series>
        void cache_series_terms(const Series &s,
        std::valarray<typename series_mult_rep<typename Series::term_type::cf_type>::type> &cfs,
        std::valarray<typename series_mult_rep<typename Series::term_type::key_type>::type> &keys)
      {
        typedef typename Series::const_sorted_iterator const_sorted_iterator;
        const const_sorted_iterator it_f = s.template nth_index<0>().end();
        size_t i=0;
        for (const_sorted_iterator it = s.template nth_index<0>().begin(); it != it_f; ++it)
        {
          series_mult_rep<typename Series::term_type::cf_type>::assign(cfs[i],it->m_cf); 
          series_mult_rep<typename Series::term_type::key_type>::assign(keys[i],it->m_key);
          ++i;
        }
      }
      // Perform plain multiplication.
      void plain_multiplication()
      {
        typedef typename term_type::multiplication_result mult_res;
        const size_t size1 = m_cfs1.size(), size2 = m_cfs2.size();
        p_assert(size1 == m_keys1.size() and size2 == m_keys2.size());
        mult_res res;
        for (size_t i = 0; i < size1; ++i)
        {
          for (size_t j = 0; j < size2; ++j)
          {
            if (m_trunc.skip(
              series_mult_rep<cf_type1>::get(m_cfs1[i]),
              series_mult_rep<key_type>::get(m_keys1[i]),
              series_mult_rep<cf_type2>::get(m_cfs2[j]),
              series_mult_rep<key_type>::get(m_keys2[j]),
              *this))
            {
              break;
            }
            term_type::multiply(
              series_mult_rep<cf_type1>::get(m_cfs1[i]),
              series_mult_rep<key_type>::get(m_keys1[i]),
              series_mult_rep<cf_type2>::get(m_cfs2[j]),
              series_mult_rep<key_type>::get(m_keys2[j]),
              res,m_args_tuple);
            insert_multiplication_result<mult_res>::run(res,*this);
          }
        }
      }
      // After the multiplication has been performed and the result stored in the temporary hash table,
      // fetch the terms from there and put them into retval.
      void plain_insert_result_into_retval()
      {
        typedef typename mult_set::const_iterator hash_iterator;
        typedef typename Series1::sorted_iterator sorted_iterator;
        term_type term;
        sorted_iterator it_hint = m_retval.template nth_index<0>().end();
        const hash_iterator it_f = m_set.end();
        for (hash_iterator it = m_set.begin(); it != it_f; ++it)
        {
          term.m_cf = it->m_cf;
          term.m_key = it->m_key;
          it_hint = m_retval.template insert<false,true>(term,m_args_tuple,it_hint);
        }
      }
      /// Adjust the sizes of the vectors representing internally the input series.
      void adjust_input_sizes()
      {
        m_cfs1.resize(m_size1);
        m_cfs2.resize(m_size2);
        m_keys1.resize(m_size1);
        m_keys2.resize(m_size2);
      }
    public:
      // These needs to be public since they will be accessed by the truncators.
      // References to the series.
      const Series1           &m_s1;
      const Series2           &m_s2;
      // Reference to the arguments tuple.
      const ArgsTuple         &m_args_tuple;
    protected:
      // Sizes of the series.
      const size_t            m_size1;
      const size_t            m_size2;
      Series1                 &m_retval;
      // Vectors of input coefficients converted for representation during series multiplication.
      std::valarray<sm_cf1>   m_cfs1;
      std::valarray<sm_cf2>   m_cfs2;
      // Vectors of input keys converted for representation during series multiplication.
      std::valarray<sm_key>   m_keys1;
      std::valarray<sm_key>   m_keys2;
      // Container to store the result of the multiplications.
      mult_set                m_set;
      // Truncator. This must be the last one defined because it will take *this
      // as parameter for construction.
      truncator_type          m_trunc;
  };
}

#endif
