/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_PSERIES_GL_REP_H
#define PIRANHA_PSERIES_GL_REP_H

#include <algorithm> // For sorting of vectors.
#include <boost/integer_traits.hpp>
// We need gmp to do arithmetics on ranges.
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>

#include "common_typedefs.h" // For max_fast_int.

namespace piranha
{
/// Generalized lexicographic representation for two pseries.
/**
 * To be used in pseries multiplication.
 */
  template <class Ps1, class Ps2>
    class pseries_gl_rep
  {
      typedef typename Ps1::trig_type::value_type mult_type;
      typedef std::valarray<std::pair<mult_type,mult_type> > e_minmax_type;
      typedef boost::integer_traits<max_fast_int> traits;
      typedef typename Ps1::ancestor::const_iterator iterator1;
      typedef typename Ps2::ancestor::const_iterator iterator2;
      typedef typename Ps1::ancestor::cf_type cf_type1;
      typedef typename Ps2::ancestor::cf_type cf_type2;
/// Coded term structure.
      template <class T>
        struct coded_term_type
      {
        mutable T     cf;
        max_fast_int  code;
        bool          flavour;
      };
/// Compact coded term structure.
/**
 * With respect to coded_term_type it does not contain flavour because this will be used in series multiplication when we
 * will have two different containers for sines and cosines.
 */
      template <class T>
        struct compact_coded_term_type
      {
        mutable T     cf;
        max_fast_int  code;
      };
    public:
      typedef coded_term_type<cf_type1> ct_type1;
      typedef coded_term_type<cf_type2> ct_type2;
      typedef compact_coded_term_type<cf_type1> cct_type1;
      typedef compact_coded_term_type<cf_type2> cct_type2;
      typedef std::valarray<ct_type1> coded_series_type1;
      typedef std::valarray<ct_type2> coded_series_type2;
      struct cct_hasher
      {
        size_t operator()(const cct_type1 &cct) const
        {
          return max_fast_int_hash(cct.code);
        }
      };
      struct cct_equal_to
      {
        bool operator()(const cct_type1 &cct1, const cct_type1 &cct2) const
        {
          return (cct1.code == cct2.code);
        }
      };
      pseries_gl_rep(const Ps1 &a, const Ps2 &b):p1(a),p2(b),
        twidth(a.trig_width()),e_minmax(twidth),viable(false),coding_vector(twidth+1)
      {
        find_minmax();
        check_viable();
// If representation is viable, let's code the series.
        if (is_viable())
        {
          code_series();
        }
      }
/// Is this representation useful?
/**
 * It is useful when code range stays inside the range of piranha::max_fast_int.
 */
      const bool &is_viable() const
      {
        return viable;
      }
/// Access first coded series.
      const coded_series_type1 &g1() const
      {
        return cs1;
      }
/// Access second coded series.
      const coded_series_type2 &g2() const
      {
        return cs2;
      }
/// Decode multiindex into external array container.
/**
 * v must support random access through operator[].
 */
      template <class Array>
        void decode_multiindex(const max_fast_int &n, Array &v) const
      {
        p_assert(twidth == v.size());
        const max_fast_int tmp = n - h_min;
        for (trig_size_t i=0;i<twidth;++i)
        {
          v[i]=((tmp%coding_vector[i+1])/(coding_vector[i])+e_minmax[i].first);
        }
      }
      const max_fast_int &g_h_min() const
      {
        return h_min;
      }
      const max_fast_int &g_h_max() const
      {
        return h_max;
      }
    private:
// Make this private to make sure we do not call default ctor.
      pseries_gl_rep()
      {}
/// Find minimum and maximum values for multipliers after multiplication.
      void find_minmax()
      {
        e_minmax_type limits1(twidth), limits2(twidth);
        const iterator1 it1_f = p1.end();
        const iterator2 it2_f = p2.end();
        iterator1 it1 = p1.begin();
        iterator2 it2 = p2.begin();
// Fill first minmax vector. This works because at this point we are sure both series have
// at least one term.
        p_assert(p1.length() >= 1 && p2.length() >= 1);
        for (trig_size_t i=0;i<twidth;++i)
        {
          limits1[i].first=limits1[i].second=it1->g_trig()->at(i);
          limits2[i].first=limits2[i].second=it2->g_trig()->at(i);
        }
        mult_type tmp;
        for (;it1!=it1_f;++it1)
        {
          for (trig_size_t j=0;j<twidth;++j)
          {
            tmp = it1->g_trig()->at(j);
            if (tmp < limits1[j].first)
              limits1[j].first = tmp;
            if (tmp > limits1[j].second)
              limits1[j].second = tmp;
          }
        }
        for (;it2!=it2_f;++it2)
        {
          for (trig_size_t j=0;j<twidth;++j)
          {
            tmp = it2->g_trig()->at(j);
            if (tmp < limits2[j].first)
              limits2[j].first = tmp;
            if (tmp > limits2[j].second)
              limits2[j].second = tmp;
          }
        }
        std::valarray<mult_type> tmp_vec(4);
        for (trig_size_t j=0;j<twidth;++j)
        {
          tmp_vec[0]=limits1[j].second+limits2[j].second;
          tmp_vec[1]=limits1[j].first+limits2[j].first;
          tmp_vec[2]=limits1[j].second-limits2[j].first;
          tmp_vec[3]=limits1[j].first-limits2[j].second;
          std::sort(&tmp_vec[0], &tmp_vec[0] + 4);
          e_minmax[j].first=tmp_vec[0];
          e_minmax[j].second=tmp_vec[3];
        }
//         for (trig_size_t j=0;j<twidth;++j)
//         {
//           std::cout << (int)e_minmax[j].first << ',' << (int)e_minmax[j].second << '\t';
//         }
//         std::cout << std::endl;
      }
/// Check whether representation is usable or not.
/**
 * Sets the viable bool flag.
 */
      void check_viable()
      {
// We must do the computations with arbitrary integers to avoid exceeding range.
        mpz_class hmin=0, hmax=0, ck=1;
        for (trig_size_t i=0;i<twidth;++i)
        {
          hmin+=ck*e_minmax[i].first;
          hmax+=ck*e_minmax[i].second;
// Assign also the coding vector, so we avoid doing it later.
          coding_vector[i]=ck.get_si();
          ck*=(e_minmax[i].second-e_minmax[i].first+1);
        }
// We want to fill on extra slot of the coding vector (wrt to "nominal" width twidth).
// This is handy for decodification.
        coding_vector[twidth]=ck.get_si();
        if (ck > traits::min() && ck < traits::max())
        {
          viable = true;
          h_min = hmin.get_si();
          h_max = hmax.get_si();
// Debug
//           std::cout << "Coding vector: ";
//           for (trig_size_t i=0;i<twidth;++i)
//           {
//             std::cout << coding_vector[i] << '\t';
//           }
//           std::cout << "+\t" << coding_vector[twidth] << '\n';
        }
// This is debug and not really needed.
        std::cout << "h: " << h_min << ',' << h_max << '\n';
      }
/// Code the series.
      void code_series()
      {
        const iterator1 it1_f = p1.end();
        const iterator2 it2_f = p2.end();
        iterator1 it1 = p1.begin();
        iterator2 it2 = p2.begin();
        const size_t l1 = p1.length(), l2 = p2.length();
        cs1.resize(l1);
        cs2.resize(l2);
        size_t i;
        for (i=0;i<l1;++i)
        {
          cs1[i].cf=*it1->g_cf();
          code_multiindex(*it1->g_trig(),cs1[i].code);
          cs1[i].flavour=it1->g_flavour();
          ++it1;
        }
        for (i=0;i<l2;++i)
        {
          cs2[i].cf=*it2->g_cf();
          code_multiindex(*it2->g_trig(),cs2[i].code);
          cs2[i].flavour=it2->g_flavour();
          ++it2;
        }
      }
/// Code a single multiindex.
      template <class T>
        void code_multiindex(const T &m,max_fast_int &n)
      {
        n=0;
        for (trig_size_t i=0;i<twidth;++i)
        {
          n+=(coding_vector[i]*m.at(i));
        }
      }
    private:
      const Ps1                               &p1;
      const Ps2                               &p2;
      const trig_size_t                       twidth;
      e_minmax_type                           e_minmax;
      bool                                    viable;
      std::valarray<max_fast_int>             coding_vector;
      coded_series_type1                      cs1;
      coded_series_type2                      cs2;
      max_fast_int                            h_min;
      max_fast_int                            h_max;
      static const boost::hash<max_fast_int>  max_fast_int_hash;
  };

  template <class Ps1, class Ps2>
    const boost::hash<max_fast_int> pseries_gl_rep<Ps1,Ps2>::max_fast_int_hash = boost::hash<max_fast_int>();
}

#endif
