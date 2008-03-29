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

#ifndef PIRANHA_FOURIER_SERIES_MULTIPLIER_H
#define PIRANHA_FOURIER_SERIES_MULTIPLIER_H

#include <boost/algorithm/minmax_element.hpp> // To calculate limits of multiplication.
#include <boost/integer_traits.hpp> // For integer limits.
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <valarray>
#include <vector>

#include "../base_classes/plain_series_multiplier.h"
#include "../integer_typedefs.h"

namespace piranha
{
  /// Series multiplier specifically tuned for Fourier series.
  /**
   * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
   * like piranha::plain_series_multiplier.
   */
  template <class Series1, class Series2, class ArgsTuple, template <class, class, class> class Truncator>
    class fourier_series_multiplier:public plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator>
  {
      typedef plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator> ancestor;
      typedef typename ancestor::truncator_type truncator_type;
      typedef typename ancestor::term_type term_type;
      typedef typename ancestor::cf_type1 cf_type1;
      typedef typename ancestor::cf_type2 cf_type2;
      typedef typename ancestor::key_type key_type;
      typedef typename ancestor::sm_cf1 sm_cf1;
      typedef typename ancestor::sm_cf2 sm_cf2;
      typedef typename ancestor::sm_key sm_key;
      typedef typename ancestor::mult_set mult_set;
      typedef boost::integer_traits<max_fast_int> traits;
    public:
      fourier_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        ancestor::plain_series_multiplier(s1,s2,retval,args_tuple),
        m_size(args_tuple.template get<key_type::position>().size()),
        m_cr_is_viable(false),
        m_min_max1(m_size),
        m_min_max2(m_size),
        m_res_min_max(m_size),
        // Coding vector is larger to accomodate extra element at the end.
        m_coding_vector(m_size+1)
      {}
      /// Perform multiplication and place the result into m_retval.
      void perform_multiplication()
      {
        find_input_min_max();
        calculate_result_min_max();
        determine_viability();
        if (m_cr_is_viable)
        {

        }
        else
        {
          ancestor::adjust_input_sizes();
          ancestor::cache_series_terms(ancestor::m_s1,ancestor::m_cfs1,ancestor::m_keys1);
          ancestor::cache_series_terms(ancestor::m_s2,ancestor::m_cfs2,ancestor::m_keys2);
          ancestor::plain_multiplication();
          ancestor::plain_insert_result_into_retval();
        }
      }
    private:
      void find_input_min_max()
      {
        typedef typename Series1::const_sorted_iterator iterator1;
        typedef typename Series2::const_sorted_iterator iterator2;
        const iterator1 it_f1 = ancestor::m_s1.template nth_index<0>().end();
        const iterator2 it_f2 = ancestor::m_s2.template nth_index<0>().end();
        iterator1 it1 = ancestor::m_s1.template nth_index<0>().begin();
        iterator2 it2 = ancestor::m_s2.template nth_index<0>().begin();
        // Fill first minmax vector. This works because at this point we are sure both series have
        // at least one term. Assert it, just to make sure.
        p_assert(!ancestor::m_s1.template nth_index<0>().empty() and !ancestor::m_s2.template nth_index<0>().empty());
        it1->m_key.template update_limits<true>(m_min_max1);
        it2->m_key.template update_limits<true>(m_min_max2);
        // Move to the second terms and cycle on all remaining terms.
        ++it1;
        ++it2;
        for (; it1 != it_f1; ++it1)
        {
          it1->m_key.template update_limits<false>(m_min_max1);
        }
        for (; it2 != it_f2; ++it2)
        {
          it2->m_key.template update_limits<false>(m_min_max2);
        }
std::cout << "Limits are:\n";
for (size_t i = 0; i < m_min_max1.size(); ++i)
{
  std::cout << m_min_max1[i].first << ',' << m_min_max1[i].second << '\n';
}
std::cout << "and:\n";
for (size_t i = 0; i < m_min_max2.size(); ++i)
{
  std::cout << m_min_max2[i].first << ',' << m_min_max2[i].second << '\n';
}
      }
      void calculate_result_min_max()
      {
        std::vector<mpz_class> tmp_vec(8);
        std::pair<typename std::vector<mpz_class>::const_iterator,std::vector<mpz_class>::const_iterator> min_max;
        for (size_t i=0; i < m_size; ++i)
        {
          tmp_vec[0]=mpz_class(m_min_max1[i].second)+mpz_class(m_min_max2[i].second);
          tmp_vec[1]=mpz_class(m_min_max1[i].first)+mpz_class(m_min_max2[i].first);
          tmp_vec[2]=mpz_class(m_min_max1[i].second)-mpz_class(m_min_max2[i].first);
          tmp_vec[3]=mpz_class(m_min_max1[i].first)-mpz_class(m_min_max2[i].second);
          tmp_vec[4]=mpz_class(m_min_max1[i].first);
          tmp_vec[5]=mpz_class(m_min_max2[i].first);
          tmp_vec[6]=mpz_class(m_min_max1[i].second);
          tmp_vec[7]=mpz_class(m_min_max2[i].second);
          min_max = boost::minmax_element(tmp_vec.begin(),tmp_vec.end());
          m_res_min_max[i].first=*(min_max.first);
          m_res_min_max[i].second=*(min_max.second);
        }
std::cout << "Mult limits are:\n";
for (size_t i = 0; i < m_res_min_max.size(); ++i)
{
  std::cout << m_res_min_max[i].first << ',' << m_res_min_max[i].second << '\n';
}
      }
      void determine_viability()
      {
        // We must do the computations with arbitrary integers to avoid exceeding range.
        mpz_class hmin(0), hmax(0), ck(1);
        for (size_t i=0; i < m_size; ++i)
        {
          hmin+=ck*m_res_min_max[i].first;
          hmax+=ck*m_res_min_max[i].second;
          // Assign also the coding vector, so we avoid doing it later.
          m_coding_vector[i]=ck.get_si();
          ck*=(m_res_min_max[i].second-m_res_min_max[i].first+1);
        }
        // We want to fill on extra slot of the coding vector (wrt to the nominal size,
        // corresponding to the arguments number for the key). This is handy for decodification.
        m_coding_vector[m_size]=ck.get_si();
        if (ck > traits::min() && ck < traits::max())
        {
          m_cr_is_viable = true;
          m_h_min = hmin.get_si();
          m_h_max = hmax.get_si();
std::cout << "Coding vector: ";
for (size_t i=0; i < m_size; ++i)
{
  std::cout << m_coding_vector[i] << '\t';
}
std::cout << "+\t" << m_coding_vector[m_size] << '\n';
        }
      }
      /// Code the series.
      void code_series()
      {
//         const iterator1 it1_f = p1.end();
//         const iterator2 it2_f = p2.end();
//         iterator1 it1 = p1.begin();
//         iterator2 it2 = p2.begin();
//         const size_t l1 = p1.length(), l2 = p2.length();
//         cs1.resize(l1);
//         cs2.resize(l2);
//         size_t i;
//         for (i=0;i<l1;++i)
//         {
//           cs1[i].cf()=it1->cf();
//           code_multiindex(it1->trig(),cs1[i].code);
//           cs1[i].flavour=it1->trig().flavour();
//           ++it1;
//         }
//         for (i=0;i<l2;++i)
//         {
//           cs2[i].cf()=it2->cf();
//           code_multiindex(it2->trig(),cs2[i].code);
//           cs2[i].flavour=it2->trig().flavour();
//           ++it2;
//         }
      }
    private:
      // Size of limits vectors (corresponding to the size of arguments vector of the key).
      const size_t                                          m_size;
      // Is coded representation viable?
      bool                                                  m_cr_is_viable;
      // Vectors of minimum and maximum value pairs for the series being multiplied.
      std::valarray<std::pair<max_fast_int,max_fast_int> >  m_min_max1;
      std::valarray<std::pair<max_fast_int,max_fast_int> >  m_min_max2;
      // Vector of minimum and maximum value pairs for the resulting series.
      // GMP is used to avoid trespassing the range limits of max_fast_int.
      std::valarray<std::pair<mpz_class,mpz_class> >        m_res_min_max;
      // Coding vector.
      std::valarray<max_fast_int>                           m_coding_vector;
      max_fast_int                                          m_h_min;
      max_fast_int                                          m_h_max;
  };
}

#endif
