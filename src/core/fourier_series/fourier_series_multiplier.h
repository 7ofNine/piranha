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
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <valarray>
#include <vector>

#include "../base_classes/plain_series_multiplier.h"
#include "../integer_typedefs.h"
#include "../memory.h"

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
      typedef typename ancestor::mult_set mult_set;
      typedef boost::integer_traits<max_fast_int> traits;
      typedef typename Series1::const_sorted_iterator iterator1;
      typedef typename Series2::const_sorted_iterator iterator2;
      typedef std::pair<cf_type1,bool> vc_res_type;
    public:
      fourier_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        ancestor::plain_series_multiplier(s1,s2,retval,args_tuple),
        m_size(args_tuple.template get<key_type::position>().size()),
        m_cr_is_viable(false),
        m_min_max1(m_size),
        m_min_max2(m_size),
        m_res_min_max(m_size),
        // Coding vector is larger to accomodate extra element at the end.
        m_coding_vector(m_size+1),
        m_vc_res(0)
      {}
      ~fourier_series_multiplier()
      {
        piranha_free(m_vc_res);
      }
      /// Perform multiplication and place the result into m_retval.
      void perform_multiplication()
      {
        find_input_min_max();
        calculate_result_min_max();
        determine_viability();
        if (m_cr_is_viable)
        {
          code_series();
          if (!perform_vector_coded_multiplication())
          {
            
          }
        }
        else
        {
          ancestor::perform_plain_multiplication();
        }
      }
    private:
      void find_input_min_max()
      {
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
          // Downcast minimum and maximum result values to fast integers.
          m_fast_res_min_max.resize(m_size);
          for (size_t i = 0; i < m_size; ++i)
          {
            m_fast_res_min_max[i].first = m_res_min_max[i].first.get_si();
            m_fast_res_min_max[i].second = m_res_min_max[i].second.get_si();
          }
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
        const iterator1 it1_f = ancestor::m_s1.template nth_index<0>().end();
        const iterator2 it2_f = ancestor::m_s2.template nth_index<0>().end();
        iterator1 it1 = ancestor::m_s1.template nth_index<0>().begin();
        iterator2 it2 = ancestor::m_s2.template nth_index<0>().begin();
        const size_t size1 = ancestor::m_s1.template nth_index<0>().size(), size2 = ancestor::m_s2.template nth_index<0>().size();
        // Make space in the coefficients and coded keys vectors.
        ancestor::m_cfs1.resize(size1);
        ancestor::m_cfs2.resize(size2);
        m_ckeys1.resize(size1);
        m_ckeys2.resize(size2);
        size_t i;
        for (i = 0; i < size1; ++i)
        {
          series_mult_rep<cf_type1>::assign(ancestor::m_cfs1[i],it1->m_cf);
          it1->m_key.code(m_coding_vector,m_ckeys1[i].first,ancestor::m_args_tuple);
          ++it1;
        }
        for (i = 0; i < size2; ++i)
        {
          series_mult_rep<cf_type2>::assign(ancestor::m_cfs2[i],it2->m_cf);
          it2->m_key.code(m_coding_vector,m_ckeys2[i].first,ancestor::m_args_tuple);
          ++it2;
        }
      }
      bool perform_vector_coded_multiplication()
      {
        // Try to allocate the space for vector coded multiplication. We need twice the number of possible
        // codes to accomodate sines and cosines results.
        // The +1 is needed because we need the number of possible codes between min and max, e.g.:
        // m_h_min = 1, m_h_max = 2 --> n of codes = 2.
        const size_t n_codes = (size_t)(m_h_max - m_h_min + 1);
        try
        {
          m_vc_res = (vc_res_type *)piranha_malloc(sizeof(vc_res_type)*(n_codes << 1));
          // Reset memory area. Use positional new so that if cf is a class with non-trivial ctors,
          // we are sure it will be initialized properly. We want to make sure the coefficients are initialized
          // to zero in order to accumulate Poisson terms during multiplication.
          for (size_t i = 0; i < (n_codes << 1); ++i)
          {
            ::new(&((m_vc_res+i)->first)) cf_type1(0,ancestor::m_args_tuple);
            (m_vc_res+i)->second = false;
          }
        }
        catch(std::bad_alloc)
        {
          return false;
        }
        // Define the base pointers for storing the results of multiplication.
        // Please note that even if here it seems like we are going to write outside allocated memory,
        // the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
        vc_res_type *vc_res_cos =  m_vc_res - m_h_min, *vc_res_sin = m_vc_res + n_codes - m_h_min;
        // Perform multiplication.
        const size_t size1 = ancestor::m_s1.template nth_index<0>().size(), size2 = ancestor::m_s2.template nth_index<0>().size();
        cf_type1 tmp_cf;
        for (size_t i = 0; i < size1; ++i)
        {
          for (size_t j = 0; j < size2; ++j)
          {
            tmp_cf = series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]);
            tmp_cf.mult_by(series_mult_rep<cf_type1>::get(ancestor::m_cfs2[j]),ancestor::m_args_tuple);
            tmp_cf.divide_by(2,ancestor::m_args_tuple);
            const max_fast_int index_plus = m_ckeys1[i].first + m_ckeys2[j].first,
              index_minus = m_ckeys1[i].first - m_ckeys2[j].first;
            switch (m_ckeys1[i].second == m_ckeys2[j].second)
            {
              case true:
                switch (m_ckeys1[i].second)
                {
                  case true:
                    vc_res_cos[index_minus].first.add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_cos[index_minus].second = true;
                    vc_res_cos[index_plus].first.add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_cos[index_plus].second = true;
                    break;
                  case false:
                    vc_res_cos[index_minus].first.add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_cos[index_minus].second = true;
                    vc_res_cos[index_plus].first.subtract(tmp_cf,ancestor::m_args_tuple);
                    vc_res_cos[index_plus].second = true;
                }
                break;
              case false:
                switch (m_ckeys1[i].second)
                {
                  case true:
                    vc_res_sin[index_minus].first.subtract(tmp_cf,ancestor::m_args_tuple);
                    vc_res_sin[index_minus].second = true;
                    vc_res_sin[index_plus].first.add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_sin[index_plus].second = true;
                    break;
                  case false:
                    vc_res_sin[index_minus].first.add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_sin[index_minus].second = true;
                    vc_res_sin[index_plus].first.add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_sin[index_plus].second = true;
                }
            }
          }
        }
        // Decode and insert the results into return value.
        term_type tmp_term;
        iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
        for (max_fast_int i = m_h_min; i <= m_h_max; ++i)
        {
          switch (unlikely(vc_res_cos[i].second))
          {
            case true:
              tmp_term.m_cf = vc_res_cos[i].first;
              tmp_term.m_key.decode(i,m_coding_vector,m_h_min,m_fast_res_min_max,ancestor::m_args_tuple);
              tmp_term.m_key.flavour() = true;
              it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
              break;
            case false:
              ;
          }
        }
        for (max_fast_int i = m_h_min; i <= m_h_max; ++i)
        {
          switch (unlikely(vc_res_sin[i].second))
          {
            case true:
              tmp_term.m_cf = vc_res_sin[i].first;
              tmp_term.m_key.decode(i,m_coding_vector,m_h_min,m_fast_res_min_max,ancestor::m_args_tuple);
              tmp_term.m_key.flavour() = false;
              it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
              break;
            case false:
              ;
          }
        }
        // Call dtors for the coefficients in the allocated space.
        // This is necessary for non-trivial coefficients.
        for (size_t i = 0; i < (n_codes << 1); ++i)
        {
          (m_vc_res+i)->first.~cf_type1();
        }
        return true;
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
      // Version of the above downcast to fast integer type.
      std::valarray<std::pair<max_fast_int,max_fast_int> >  m_fast_res_min_max;
      // Coding vector.
      std::valarray<max_fast_int>                           m_coding_vector;
      max_fast_int                                          m_h_min;
      max_fast_int                                          m_h_max;
      // Coded keys.
      std::valarray<typename key_type::coded_type>          m_ckeys1;
      std::valarray<typename key_type::coded_type>          m_ckeys2;
      // Result of vector coded multiplication.
      vc_res_type                                           *m_vc_res;
  };
}

#endif
