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

#ifndef PIRANHA_POLYNOMIAL_MULTIPLIER_H
#define PIRANHA_POLYNOMIAL_MULTIPLIER_H

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
  /// Series multiplier specifically tuned for Polynomials.
  /**
   * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
   * like piranha::plain_series_multiplier.
   */
  template <class Series1, class Series2, class ArgsTuple, template <class, class, class> class Truncator>
    class polynomial_multiplier:
    public plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator>,
    public coded_series_multiplier<polynomial_multiplier<Series1,Series2,ArgsTuple,Truncator> >
  {
      typedef plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator> ancestor;
      typedef coded_series_multiplier<polynomial_multiplier<Series1,Series2,ArgsTuple,Truncator> > coded_ancestor;
      friend class coded_series_multiplier<polynomial_multiplier<Series1,Series2,ArgsTuple,Truncator> >;
      // Some of these typedefs are used also in the coded ancestor.
      typedef typename ancestor::truncator_type truncator_type;
      typedef typename ancestor::term_type term_type;
      typedef typename ancestor::cf_type1 cf_type1;
      typedef typename ancestor::cf_type2 cf_type2;
      typedef typename ancestor::key_type key_type;
      typedef typename ancestor::mult_set mult_set;
      typedef typename Series1::const_sorted_iterator iterator1;
      typedef typename Series2::const_sorted_iterator iterator2;
    public:
      polynomial_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
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
        std::vector<mpz_class> tmp_vec(6);
        std::pair<typename std::vector<mpz_class>::const_iterator,std::vector<mpz_class>::const_iterator> min_max;
        for (size_t i=0; i < coded_ancestor::m_size; ++i)
        {
          tmp_vec[0]=mpz_class(coded_ancestor::m_min_max1[i].second)+mpz_class(coded_ancestor::m_min_max2[i].second);
          tmp_vec[1]=mpz_class(coded_ancestor::m_min_max1[i].first)+mpz_class(coded_ancestor::m_min_max2[i].first);
          tmp_vec[2]=mpz_class(coded_ancestor::m_min_max1[i].first);
          tmp_vec[3]=mpz_class(coded_ancestor::m_min_max2[i].first);
          tmp_vec[4]=mpz_class(coded_ancestor::m_min_max1[i].second);
          tmp_vec[5]=mpz_class(coded_ancestor::m_min_max2[i].second);
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
      bool perform_vector_coded_multiplication()
      {
        cf_type1 *p_vc_res(0);
        // Try to allocate the space for vector coded multiplication.
        // The +1 is needed because we need the number of possible codes between min and max, e.g.:
        // coded_ancestor::m_h_min = 1, coded_ancestor::m_h_max = 2 --> n of codes = 2.
        const size_t n_codes = (size_t)(coded_ancestor::m_h_max - coded_ancestor::m_h_min + 1);
        try
        {
          p_vc_res = (cf_type1 *)piranha_malloc(sizeof(cf_type1)*n_codes);
          // Reset memory area. Use positional new so that if cf is a class with non-trivial ctors,
          // we are sure it will be initialized properly. We want to make sure the coefficients are initialized
          // to zero in order to accumulate monomials during multiplication.
          for (size_t i = 0; i < n_codes; ++i)
          {
            ::new(p_vc_res+i) cf_type1(0,ancestor::m_args_tuple);
          }
        }
        catch(std::bad_alloc)
        {
          piranha_free(p_vc_res);
          return false;
        }
        // Define the base pointers for storing the results of multiplication.
        // Please note that even if here it seems like we are going to write outside allocated memory,
        // the indices from the analysis of the coded series will prevent out-of-boundaries reads/writes.
        cf_type1 *vc_res =  p_vc_res - coded_ancestor::m_h_min;
        // Perform multiplication.
        for (size_t i = 0; i < ancestor::m_size1; ++i)
        {
          for (size_t j = 0; j < ancestor::m_size2; ++j)
          {
            if (ancestor::m_trunc.skip(
              series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]),
              series_mult_rep<key_type>::get(ancestor::m_keys1[i]),
              series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),
              series_mult_rep<key_type>::get(ancestor::m_keys2[j])))
            {
              break;
            }
            // TODO: Does it make sense here to define a method for coefficients like:
            // mult_by_and_insert_into<bool Sign>(cf2,retval,m_args_tuple)
            // so that we can avoid copying stuff around here and elsewhere?
            cf_type1 tmp_cf(series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]));
            tmp_cf.mult_by(series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),ancestor::m_args_tuple);
            vc_res[coded_ancestor::m_ckeys1[i] + coded_ancestor::m_ckeys2[j]].add(tmp_cf,ancestor::m_args_tuple);
          }
        }
        // Decode and insert the results into return value.
        term_type tmp_term;
        iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
        for (max_fast_int i = coded_ancestor::m_h_min; i <= coded_ancestor::m_h_max; ++i)
        {
          // Take a shortcut and check for ignorability of the coefficient here.
          // This way we avoid decodification, and all the series term insertion yadda-yadda.
          switch (likely(vc_res[i].is_ignorable(ancestor::m_args_tuple)))
          {
            case true:
              break;
            case false:
              tmp_term.m_cf = vc_res[i];
              tmp_term.m_key.decode(i,coded_ancestor::m_coding_vector,coded_ancestor::m_h_min,
                coded_ancestor::m_fast_res_min_max,ancestor::m_args_tuple);
              it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
          }
        }
        // Call dtors for the coefficients in the allocated space.
        // This is necessary for non-trivial coefficients.
        for (size_t i = 0; i < n_codes; ++i)
        {
          (p_vc_res+i)->~cf_type1();
        }
        // Free the allocated space.
        piranha_free(p_vc_res);
        return true;
      }
      void perform_hash_coded_multiplication()
      {
        typedef typename coded_ancestor::template generic_cterm<cf_type1> cterm;
        typedef typename coded_ancestor::template generic_cmult_set<cf_type1>::type cmult_set;
        typedef typename cmult_set::iterator c_iterator;
        c_iterator it;
        cmult_set cms;
        // Set max load factors.
        cms.max_load_factor(settings_manager::get_load_factor());
        for (size_t i = 0; i < ancestor::m_size1; ++i)
        {
          for (size_t j = 0; j < ancestor::m_size2; ++j)
          {
            if (ancestor::m_trunc.skip(
              series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]),
              series_mult_rep<key_type>::get(ancestor::m_keys1[i]),
              series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),
              series_mult_rep<key_type>::get(ancestor::m_keys2[j])))
            {
              break;
            }
            // TODO: can this be optimized to be run in just one pass (for example, for GMP using expression
            // templates)?
            cterm tmp_term(series_mult_rep<cf_type1>::get(ancestor::m_cfs1[i]),coded_ancestor::m_ckeys1[i]);
            // Calculate the coefficient and key.
            tmp_term.m_cf.mult_by(series_mult_rep<cf_type2>::get(ancestor::m_cfs2[j]),ancestor::m_args_tuple);
            tmp_term.m_ckey += coded_ancestor::m_ckeys2[j];
            // Now find the term and handle it.
            it = cms.find(tmp_term.m_ckey);
            switch (it == cms.end())
            {
              case true:
                cms.insert(tmp_term);
                break;
              case false:
                it->m_cf.add(tmp_term.m_cf,ancestor::m_args_tuple);
            }
          }
        }
        // Decode and insert into retval.
        term_type tmp_term;
        iterator1 it_hint = ancestor::m_retval.template nth_index<0>().end();
        const c_iterator c_it_f = cms.end();
        for (c_iterator c_it = cms.begin(); c_it != c_it_f; ++c_it)
        {
          tmp_term.m_cf = c_it->m_cf;
          tmp_term.m_key.decode(c_it->m_ckey,coded_ancestor::m_coding_vector,coded_ancestor::m_h_min,
            coded_ancestor::m_fast_res_min_max,ancestor::m_args_tuple);
          it_hint = ancestor::m_retval.insert(tmp_term,ancestor::m_args_tuple,it_hint);
        }
      }
  };
}

#endif
