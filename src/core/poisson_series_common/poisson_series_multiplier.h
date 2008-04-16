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

#ifndef PIRANHA_POISSON_SERIES_MULTIPLIER_H
#define PIRANHA_POISSON_SERIES_MULTIPLIER_H

#include <boost/algorithm/minmax_element.hpp> // To calculate limits of multiplication.
#include <exception>
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <valarray>
#include <vector>

#include "../base_classes/plain_series_multiplier.h"
#include "../base_classes/coded_series_multiplier.h"
#include "../integer_typedefs.h"
#include "../memory.h"
#include "../settings_manager.h" // For hash load factor.

namespace piranha
{
  /// Series multiplier specifically tuned for Poisson series.
  /**
   * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
   * like piranha::plain_series_multiplier.
   */
  template <class Series1, class Series2, class ArgsTuple, template <class> class Truncator>
    class poisson_series_multiplier:
    public plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator>,
    public coded_series_multiplier<poisson_series_multiplier<Series1,Series2,ArgsTuple,Truncator> >
  {
      typedef plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator> ancestor;
      typedef coded_series_multiplier<poisson_series_multiplier<Series1,Series2,ArgsTuple,Truncator> > coded_ancestor;
      friend class coded_series_multiplier<poisson_series_multiplier<Series1,Series2,ArgsTuple,Truncator> >;
      typedef typename Series1::const_sorted_iterator const_iterator1;
      typedef typename Series2::const_sorted_iterator const_iterator2;
      typedef typename Series1::sorted_iterator iterator1;
      typedef typename Series2::sorted_iterator iterator2;
    public:
      typedef typename ancestor::truncator_type truncator_type;
      typedef typename ancestor::term_type term_type;
      typedef typename ancestor::cf_type1 cf_type1;
      typedef typename ancestor::cf_type2 cf_type2;
      typedef typename ancestor::key_type key_type;
      poisson_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        ancestor::plain_series_multiplier(s1,s2,retval,args_tuple)
      {}
      /// Perform multiplication and place the result into m_retval.
      void perform_multiplication()
      {
        coded_ancestor::find_input_min_max();
        calculate_result_min_max();
        coded_ancestor::determine_viability();
        if (coded_ancestor::m_cr_is_viable)
        {
          coded_ancestor::store_coefficients_code_keys();
          store_flavours();
          if (!perform_vector_coded_multiplication())
          {
            std::cout << "Going for hashed!\n";
            perform_hash_coded_multiplication();
          }
        }
        else
        {
          ancestor::perform_plain_multiplication();
        }
      }
    private:
      void calculate_result_min_max()
      {
        // TODO: optimize here the usage of mpz classes.
        std::vector<mpz_class> tmp_vec(8);
        std::pair<typename std::vector<mpz_class>::const_iterator,std::vector<mpz_class>::const_iterator> min_max;
        for (size_t i=0; i < coded_ancestor::m_size; ++i)
        {
          tmp_vec[0]=mpz_class(coded_ancestor::m_min_max1[i].second)+mpz_class(coded_ancestor::m_min_max2[i].second);
          tmp_vec[1]=mpz_class(coded_ancestor::m_min_max1[i].first)+mpz_class(coded_ancestor::m_min_max2[i].first);
          tmp_vec[2]=mpz_class(coded_ancestor::m_min_max1[i].second)-mpz_class(coded_ancestor::m_min_max2[i].first);
          tmp_vec[3]=mpz_class(coded_ancestor::m_min_max1[i].first)-mpz_class(coded_ancestor::m_min_max2[i].second);
          tmp_vec[4]=mpz_class(coded_ancestor::m_min_max1[i].first);
          tmp_vec[5]=mpz_class(coded_ancestor::m_min_max2[i].first);
          tmp_vec[6]=mpz_class(coded_ancestor::m_min_max1[i].second);
          tmp_vec[7]=mpz_class(coded_ancestor::m_min_max2[i].second);
          min_max = boost::minmax_element(tmp_vec.begin(),tmp_vec.end());
          coded_ancestor::m_res_min_max[i].first=*(min_max.first);
          coded_ancestor::m_res_min_max[i].second=*(min_max.second);
        }
std::cout << "Mult limits are:\n";
for (size_t i = 0; i < coded_ancestor::m_res_min_max.size(); ++i)
{
  std::cout << coded_ancestor::m_res_min_max[i].first << ',' << coded_ancestor::m_res_min_max[i].second << '\n';
}
      }
      /// Store flavours of the series into own vectors.
      void store_flavours()
      {
        iterator1 it1 = ancestor::m_s1.template nth_index<0>().begin();
        iterator2 it2 = ancestor::m_s2.template nth_index<0>().begin();
        // Make space in the flavours vectors.
        m_flavours1.resize(ancestor::m_size1);
        m_flavours2.resize(ancestor::m_size2);
        size_t i;
        for (i = 0; i < ancestor::m_size1; ++i)
        {
          m_flavours1[i] = it1->m_key.flavour();
          ++it1;
        }
        for (i = 0; i < ancestor::m_size2; ++i)
        {
          m_flavours2[i] = it2->m_key.flavour();
          ++it2;
        }
      }
      bool perform_vector_coded_multiplication()
      {
        cf_type1 *p_vc_res_cos(0), *p_vc_res_sin(0);
        // Try to allocate the space for vector coded multiplication. We need two arrays of results,
        // one for cosines, one for sines.
        // The +1 is needed because we need the number of possible codes between min and max, e.g.:
        // coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
        const size_t n_codes = (size_t)(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1);
        try
        {
          p_vc_res_cos = (cf_type1 *)piranha_malloc(sizeof(cf_type1)*n_codes);
          p_vc_res_sin = (cf_type1 *)piranha_malloc(sizeof(cf_type1)*n_codes);
          // Reset memory area. Use positional new so that if cf is a class with non-trivial ctors,
          // we are sure it will be initialized properly. We want to make sure the coefficients are initialized
          // to zero in order to accumulate Poisson terms during multiplication.
          for (size_t i = 0; i < n_codes; ++i)
          {
            ::new(p_vc_res_cos+i) cf_type1(0,ancestor::m_args_tuple);
            ::new(p_vc_res_sin+i) cf_type1(0,ancestor::m_args_tuple);
          }
        }
        catch(const std::bad_alloc &)
        {
          piranha_free(p_vc_res_cos);
          piranha_free(p_vc_res_sin);
          return false;
        }
std::cout << "Going for Poisson coded\n";
        // Define the base pointers for storing the results of multiplication.
        // Please note that even if here it seems like we are going to write outside allocated memory,
        // the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
        cf_type1 *vc_res_cos =  p_vc_res_cos - coded_ancestor::m_h_min, *vc_res_sin = p_vc_res_sin - coded_ancestor::m_h_min;
        // Perform multiplication.
        for (size_t i = 0; i < ancestor::m_size1; ++i)
        {
          for (size_t j = 0; j < ancestor::m_size2; ++j)
          {
            if (ancestor::m_trunc.skip(
              series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]),
              series_mult_rep<key_type>::get(ancestor::m_keys1[i]),
              series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),
              series_mult_rep<key_type>::get(ancestor::m_keys2[j]),
              *this))
            {
              break;
            }
            // TODO: Does it make sense here to define a method for coefficients like:
            // mult_by_and_insert_into<bool Sign>(cf2,retval,m_args_tuple)
            // so that we can avoid copying stuff around here and elsewhere?
            // TODO: don't create a tmp_cf each time, create it outside the loop and assign it here.
            // For nontrivial coefficients we can save a lot of memory allocations.
            cf_type1 tmp_cf(series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]));
            tmp_cf.mult_by(series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),ancestor::m_args_tuple);
            tmp_cf.divide_by(2,ancestor::m_args_tuple);
            const max_fast_int index_plus = coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j],
              index_minus = coded_ancestor::m_ckeys1[i] - coded_ancestor::m_ckeys2[j];
            switch (m_flavours1[i] == m_flavours2[j])
            {
              case true:
                switch (m_flavours1[i])
                {
                  case true:
                    vc_res_cos[index_minus].add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_cos[index_plus].add(tmp_cf,ancestor::m_args_tuple);
                    break;
                  case false:
                    vc_res_cos[index_minus].add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_cos[index_plus].subtract(tmp_cf,ancestor::m_args_tuple);
                }
                break;
              case false:
                switch (m_flavours1[i])
                {
                  case true:
                    vc_res_sin[index_minus].subtract(tmp_cf,ancestor::m_args_tuple);
                    vc_res_sin[index_plus].add(tmp_cf,ancestor::m_args_tuple);
                    break;
                  case false:
                    vc_res_sin[index_minus].add(tmp_cf,ancestor::m_args_tuple);
                    vc_res_sin[index_plus].add(tmp_cf,ancestor::m_args_tuple);
                }
            }
          }
        }
        // Decode and insert the results into return value.
        term_type tmp_term;
        iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
        for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i)
        {
          // Take a shortcut and check for ignorability of the coefficient here.
          // This way we avoid decodification, and all the series term insertion yadda-yadda.
          switch (likely(vc_res_cos[i].is_ignorable(ancestor::m_args_tuple)))
          {
            case true:
              break;
            case false:
              tmp_term.m_cf = vc_res_cos[i];
              coded_ancestor::decode(tmp_term.m_key,i);
              tmp_term.m_key.flavour() = true;
              it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
          }
        }
        for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i)
        {
          switch (likely(vc_res_sin[i].is_ignorable(ancestor::m_args_tuple)))
          {
            case true:
              break;
            case false:
              tmp_term.m_cf = vc_res_sin[i];
              coded_ancestor::decode(tmp_term.m_key,i);
              tmp_term.m_key.flavour() = false;
              it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
          }
        }
        // Call dtors for the coefficients in the allocated space.
        // This is necessary for non-trivial coefficients.
        for (size_t i = 0; i < n_codes; ++i)
        {
          (p_vc_res_cos+i)->~cf_type1();
          (p_vc_res_sin+i)->~cf_type1();
        }
        // Free the allocated space.
        piranha_free(p_vc_res_cos);
        piranha_free(p_vc_res_sin);
        return true;
      }
      void perform_hash_coded_multiplication()
      {
        typedef typename coded_ancestor::template generic_cterm<cf_type1> cterm;
        typedef typename coded_ancestor::template generic_cmult_set<cf_type1>::type cmult_set;
        typedef typename cmult_set::iterator c_iterator;
        c_iterator it;
        cmult_set cms_cos, cms_sin;
        // Set max load factors.
        cms_cos.max_load_factor(settings_manager::load_factor());
        cms_sin.max_load_factor(settings_manager::load_factor());
        for (size_t i = 0; i < ancestor::m_size1; ++i)
        {
          for (size_t j = 0; j < ancestor::m_size2; ++j)
          {
            if (ancestor::m_trunc.skip(
              series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]),
              series_mult_rep<key_type>::get(ancestor::m_keys1[i]),
              series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),
              series_mult_rep<key_type>::get(ancestor::m_keys2[j]),
              *this))
            {
              break;
            }
            // TODO: here (and elsewhere, likely), we can avoid an extra copy by working with keys and cfs instead of terms,
            // generating only one coefficient and change its sign later if needed - after insertion.
            cterm tmp_term1(series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]),coded_ancestor::m_ckeys1[i]);
            // Handle the coefficient, with positive signs for now.
            tmp_term1.m_cf.mult_by(series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),ancestor::m_args_tuple);
            tmp_term1.m_cf.divide_by(2,ancestor::m_args_tuple);
            tmp_term1.m_ckey -= coded_ancestor::m_ckeys2[j];
            // Create the second term, using the first one's coefficient and the appropriate code.
            cterm tmp_term2(tmp_term1.m_cf,coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j]);
            // Now fix flavours and coefficient signs.
            switch (m_flavours1[i] == m_flavours2[j])
            {
              case true:
                switch (m_flavours1[i])
                {
                  case true:
                    break;
                  case false:
                    tmp_term2.m_cf.invert_sign(ancestor::m_args_tuple);
                }
                // Insert into cosine container.
                it = cms_cos.find(tmp_term1.m_ckey);
                switch (it == cms_cos.end())
                {
                  case true:
                    cms_cos.insert(tmp_term1);
                    break;
                  case false:
                    it->m_cf.add(tmp_term1.m_cf,ancestor::m_args_tuple);
                }
                it = cms_cos.find(tmp_term2.m_ckey);
                switch (it == cms_cos.end())
                {
                  case true:
                    cms_cos.insert(tmp_term2);
                    break;
                  case false:
                    it->m_cf.add(tmp_term2.m_cf,ancestor::m_args_tuple);
                }
                break;
              case false:
                switch (m_flavours1[i])
                {
                  case true:
                    tmp_term1.m_cf.invert_sign(ancestor::m_args_tuple);
                    break;
                  case false:
                    ;
                }
                // Insert into sine container.
                it = cms_sin.find(tmp_term1.m_ckey);
                switch (it == cms_sin.end())
                {
                  case true:
                    cms_sin.insert(tmp_term1);
                    break;
                  case false:
                    it->m_cf.add(tmp_term1.m_cf,ancestor::m_args_tuple);
                }
                it = cms_sin.find(tmp_term2.m_ckey);
                switch (it == cms_sin.end())
                {
                  case true:
                    cms_sin.insert(tmp_term2);
                    break;
                  case false:
                    it->m_cf.add(tmp_term2.m_cf,ancestor::m_args_tuple);
                }
            }
          }
        }
        term_type tmp_term;
        iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
        {
          const c_iterator c_it_f = cms_cos.end();
          for (c_iterator c_it = cms_cos.begin(); c_it != c_it_f; ++c_it)
          {
            tmp_term.m_cf = c_it->m_cf;
            coded_ancestor::decode(tmp_term.m_key,c_it->m_ckey);
            tmp_term.m_key.flavour() = true;
            it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
          }
        }
        {
          const c_iterator c_it_f = cms_sin.end();
          for (c_iterator c_it = cms_sin.begin(); c_it != c_it_f; ++c_it)
          {
            tmp_term.m_cf = c_it->m_cf;
            coded_ancestor::decode(tmp_term.m_key,c_it->m_ckey);
            tmp_term.m_key.flavour() = false;
            it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
          }
        }
      }
    private:
      // For Poisson series we also need flavours.
      std::valarray<bool> m_flavours1;
      std::valarray<bool> m_flavours2;
  };
}

#endif
