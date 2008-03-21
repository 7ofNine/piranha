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
      typedef typename Series1::cf_type cf1;
      typedef typename Series2::cf_type cf2;
      typedef typename Series1::key_type key;
      typedef typename convert_for_series_multiplication<cf1>::type sm_cf1;
      typedef typename convert_for_series_multiplication<cf2>::type sm_cf2;
      typedef typename convert_for_series_multiplication<key>::type sm_key;
    public:
      series_mult_rep(const Series1 &s1, const Series2 &s2):
        m_s1(s1),m_s2(s2),m_cfs1(s1.template nth_index<0>().size()),m_cfs2(s2.template nth_index<0>().size()),
        m_keys1(s1.template nth_index<0>().size()),m_keys2(s2.template nth_index<0>().size())
      {

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
