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

#ifndef PIRANHA_CODED_SERIES_MULTIPLIER_H
#define PIRANHA_CODED_SERIES_MULTIPLIER_H

#include <boost/functional/hash.hpp>
#include <boost/integer_traits.hpp> // For integer limits.
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <gmp.h>
#include <gmpxx.h>
#include <utility> // For std::pair.
#include <valarray>

#include "../p_assert.h"
#include "../settings.h" // For debug messages.

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
  /// Toolbox for coded series multiplication.
  /**
   * Intended to be inherited together with piranha::plain_series_multiplier. It adds common methods for coded
   * series multiplication.
   */
  template <class Derived>
    class coded_series_multiplier
  {
      typedef boost::integer_traits<max_fast_int> traits;
    public:
      // Decode code n into key.
      template <class Key>
        void decode(Key &key, const max_fast_int &n) const
      {
        key.decode(n,m_coding_vector,m_h_min,m_fast_res_min_max,derived_const_cast->m_args_tuple);
      }
    protected:
      // Typedefs for hash coded arithmetics.
      // This is the representation of a term for the hash coded representation.
      template <class Cf>
        struct generic_cterm
      {
        mutable Cf        m_cf;
        max_fast_int      m_ckey;
        // Templatized this way because we want to build complexes from reals.
        template <class Cf2>
          generic_cterm(const Cf2 &cf2, const max_fast_int &code):m_cf(cf2),m_ckey(code) {}
        struct hasher
        {
          size_t operator()(const max_fast_int &code) const
          {
            return boost::hash<max_fast_int>()(code);
          }
        };
        struct equal_to
        {
          bool operator()(const max_fast_int &code1, const max_fast_int &code2) const
          {
            return (code1 == code2);
          }
        };
      };
      // Template typedef for the hash set for coded multiplication.
      template <class Cf>
        struct generic_cmult_set
      {
        typedef boost::multi_index_container
        <
          generic_cterm<Cf>,
          boost::multi_index::indexed_by
          <
            boost::multi_index::hashed_unique<boost::multi_index::member<generic_cterm<Cf>,max_fast_int,&generic_cterm<Cf>::m_ckey>,
            typename generic_cterm<Cf>::hasher,typename generic_cterm<Cf>::equal_to>
          >
        > type;
      };
      coded_series_multiplier():
        m_cr_is_viable(false),
        m_size(derived_const_cast->m_args_tuple.template get<Derived::key_type::position>().size()),
        m_min_max1(m_size),m_min_max2(m_size),m_res_min_max(m_size),m_fast_res_min_max(m_size),
        // Coding vector is larger to accomodate extra element at the end.
        m_coding_vector(m_size+1)
      {}
      void find_input_min_max()
      {
        typedef typename Derived::const_iterator1 const_iterator1;
        typedef typename Derived::const_iterator2 const_iterator2;
        const const_iterator1 it_f1 = derived_const_cast->m_s1.template nth_index<0>().end();
        const const_iterator2 it_f2 = derived_const_cast->m_s2.template nth_index<0>().end();
        const_iterator1 it1 = derived_const_cast->m_s1.template nth_index<0>().begin();
        const_iterator2 it2 = derived_const_cast->m_s2.template nth_index<0>().begin();
        // Fill first minmax vector. This works because at this point we are sure both series have
        // at least one term. Assert it, just to make sure.
        p_assert(!derived_const_cast->m_s1.template nth_index<0>().empty() and
          !derived_const_cast->m_s2.template nth_index<0>().empty());
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
// std::cout << "Limits are:\n";
// for (size_t i = 0; i < m_min_max1.size(); ++i)
// {
//   std::cout << m_min_max1[i].first << ',' << m_min_max1[i].second << '\n';
// }
// std::cout << "and:\n";
// for (size_t i = 0; i < m_min_max2.size(); ++i)
// {
//   std::cout << m_min_max2[i].first << ',' << m_min_max2[i].second << '\n';
// }
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
        __PDEBUG(std::cout << "hmax-hmin=" << hmax-hmin << '\n');
        // We want to fill on extra slot of the coding vector (wrt to the nominal size,
        // corresponding to the arguments number for the key). This is handy for decodification.
        m_coding_vector[m_size]=ck.get_si();
        p_assert(ck > 0);
        // Determine viability by checking that ck and the minimum/maximum values for the codes
        // respect the fast integer boundaries.
        if (ck < traits::max() and hmin > traits::min() and hmin < traits::max() and
          hmax > traits::min() and hmax < traits::max())
        {
          m_cr_is_viable = true;
          m_h_min = hmin.get_si();
          m_h_max = hmax.get_si();
          // Downcast minimum and maximum result values to fast integers.
          for (size_t i = 0; i < m_size; ++i)
          {
            if (m_res_min_max[i].first < traits::min() or m_res_min_max[i].first > traits::max() or
              m_res_min_max[i].second < traits::min() or m_res_min_max[i].second > traits::max())
            {
              std::cout << "Warning: results of series multiplication cross " <<
              "fast integer limits. Expect errors." << std::endl;
            }
            m_fast_res_min_max[i].first = m_res_min_max[i].first.get_si();
            m_fast_res_min_max[i].second = m_res_min_max[i].second.get_si();
          }
// std::cout << "Coding vector: ";
// for (size_t i=0; i < m_size; ++i)
// {
//   std::cout << m_coding_vector[i] << '\t';
// }
// std::cout << "+\t" << m_coding_vector[m_size] << '\n';
        }
      }
      /// Store coefficients and code keys.
      void store_coefficients_code_keys()
      {
        typedef typename Derived::const_iterator1 const_iterator1;
        typedef typename Derived::const_iterator2 const_iterator2;
        typedef typename Derived::cf_type1 cf_type1;
        typedef typename Derived::cf_type2 cf_type2;
        const_iterator1 it1 = derived_const_cast->m_s1.template nth_index<0>().begin();
        const_iterator2 it2 = derived_const_cast->m_s2.template nth_index<0>().begin();
        // Make space in the coefficients and coded keys vectors.
        derived_const_cast->m_cfs1.resize(derived_const_cast->m_size1);
        derived_const_cast->m_cfs2.resize(derived_const_cast->m_size2);
        m_ckeys1.resize(derived_const_cast->m_size1);
        m_ckeys2.resize(derived_const_cast->m_size2);
        size_t i;
        for (i = 0; i < derived_const_cast->m_size1; ++i)
        {
          derived_const_cast->m_cfs1[i] = it1->m_cf;
          it1->m_key.code(m_coding_vector,m_ckeys1[i],derived_const_cast->m_args_tuple);
          ++it1;
        }
        for (i = 0; i < derived_const_cast->m_size2; ++i)
        {
          derived_const_cast->m_cfs2[i] = it2->m_cf;
          it2->m_key.code(m_coding_vector,m_ckeys2[i],derived_const_cast->m_args_tuple);
          ++it2;
        }
      }
    protected:
      // Is coded representation viable?
      bool                                                  m_cr_is_viable;
      // Size of the coding vector, min_max vectors, etc.
      const size_t                                          m_size;
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
      // Mininum and maximum values of codes.
      max_fast_int                                          m_h_min;
      max_fast_int                                          m_h_max;
      // Coded keys.
      std::valarray<max_fast_int>                           m_ckeys1;
      std::valarray<max_fast_int>                           m_ckeys2;
  };
}

#undef derived_const_cast
#undef derived_cast

#endif