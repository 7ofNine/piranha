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

#ifndef PIRANHA_SERIES_MULT_REP_H
#define PIRANHA_SERIES_MULT_REP_H

#include <valarray>

#include "../type_traits.h"

namespace piranha
{
  template <class Series1, class Series2>
    class series_mult_rep
  {
      typedef typename Series1::term_type::cf_type cf1;
      typedef typename Series2::term_type::cf_type cf2;
      typedef typename Series1::term_type::key_type key;
      typedef typename constant_copy<cf1>::type sm_cf1;
      typedef typename constant_copy<cf2>::type sm_cf2;
      typedef typename constant_copy<key>::type sm_key;
      typedef typename Series1::const_sorted_iterator const_sorted_iterator1;
      typedef typename Series2::const_sorted_iterator const_sorted_iterator2;
    public:
      series_mult_rep(const Series1 &s1, const Series2 &s2):
        m_s1(s1),m_s2(s2),m_cfs1(s1.template nth_index<0>().size()),m_cfs2(s2.template nth_index<0>().size()),
        m_keys1(s1.template nth_index<0>().size()),m_keys2(s2.template nth_index<0>().size())
      {
        cache_series_terms(s1,m_cfs1,m_keys1);
        cache_series_terms(s2,m_cfs2,m_keys2);
      }
    private:
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
    protected:
      // References to the series.
      const Series1           &m_s1;
      const Series2           &m_s2;
      // Vectors of input coefficients converted for representation during series multiplication.
      std::valarray<sm_cf1>   m_cfs1;
      std::valarray<sm_cf2>   m_cfs2;
      // Vectors of input keys converted for representation during series multiplication.
      std::valarray<sm_key>   m_keys1;
      std::valarray<sm_key>   m_keys2;
  };
}

#endif
