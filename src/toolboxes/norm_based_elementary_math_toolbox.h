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

#ifndef PIRANHA_NORM_BASED_ELEMENTARY_MATH_TOOLBOX_H
#define PIRANHA_NORM_BASED_ELEMENTARY_MATH_TOOLBOX_H

#include <algorithm> // For sorting of vectors.
#include <boost/foreach.hpp>
#include <boost/integer_traits.hpp>
#include <gmp.h>
#include <gmpxx.h>

#include "../bits/config.h" // For selection of temporary hash container for multiplication
#include "../bits/common_typedefs.h"
#include "../bits/light_term.h"

namespace piranha
{
  template <class T, class U, class MultType>
    class series_gl_rep
  {
      typedef MultType mult_type;
      typedef std::valarray<std::pair<mult_type,mult_type> > e_minmax_type;
      typedef boost::integer_traits<max_fast_int> traits;
    public:
      series_gl_rep(const T &a, const U &b, const trig_size_t &w):v_p1(a),v_p2(b),
        twidth(w),e_minmax(twidth),viable(false),coding_vector(twidth+1)
      {
        find_minmax();
        check_viable();
// If representation is viable, let's code the series.
        if (is_viable())
        {
        }
      }
      const bool &is_viable() const
      {
        return viable;
      }
    private:
// Make this private to make sure we do not call default ctor.
      series_gl_rep()
      {}
      void find_minmax()
      {
        const size_t l1 = v_p1.size(), l2 = v_p2.size();
        e_minmax_type limits1(twidth), limits2(twidth);
// Fill first minmax vector. This works because at this point we are sure both series have
// at least one term.
        p_assert(l1 >= 1 && l2 >= 1);
        for (trig_size_t i=0;i<twidth;++i)
        {
          limits1[i].first=limits1[i].second=v_p1[0]->g_trig()->at(i);
          limits2[i].first=limits2[i].second=v_p2[0]->g_trig()->at(i);
        }
        mult_type tmp;
        for (size_t i=1;i<l1;++i)
        {
          for (trig_size_t j=0;j<twidth;++j)
          {
            tmp = v_p1[i]->g_trig()->at(j);
            if (tmp < limits1[j].first)
              limits1[j].first = tmp;
            if (tmp > limits1[j].second)
              limits1[j].second = tmp;
          }
        }
        for (size_t i=1;i<l2;++i)
        {
          for (trig_size_t j=0;j<twidth;++j)
          {
            tmp = v_p2[i]->g_trig()->at(j);
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
        for (trig_size_t j=0;j<twidth;++j)
        {
          std::cout << (int)e_minmax[j].first << ',' << (int)e_minmax[j].second << '\t';
        }
        std::cout << std::endl;
      }
      void check_viable()
      {
// We must do the computations with arbitrary integers not to exceed range.
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
        if (hmin > traits::min() && hmax < traits::max())
        {
          viable = true;
// Debug
          std::cout << "Coding vector: ";
          for (trig_size_t i=0;i<twidth;++i)
          {
            std::cout << coding_vector[i] << '\t';
          }
          std::cout << "+\t" << coding_vector[twidth] << '\n';
        }
// This is debug and not really needed.
        std::cout << "h: " << hmin << ',' << hmax << '\n';
      }
    private:
      const T                     &v_p1;
      const U                     &v_p2;
      const trig_size_t           twidth;
      e_minmax_type               e_minmax;
      bool                        viable;
      std::valarray<max_fast_int> coding_vector;
  };

/// Elementary math toolbox for numerical series.
/**
 * Implements series multiplication and truncation according to norms.
 */
  template <class DerivedPs>
    class norm_based_elementary_math_toolbox
  {
    public:
      template <class Cf>
        void cf_multiplication(const Cf &cf)
      {
        typedef typename DerivedPs::ancestor::term_type term_type;
        typedef typename DerivedPs::ancestor::it_s_index it_s_index;
        DerivedPs *derived_cast=static_cast<DerivedPs *>(this);
        if (derived_cast->empty())
        {
          return;
        }
        DerivedPs tmp_ps;
        tmp_ps.merge_args(*derived_cast);
        term_type tmp_term;
        it_s_index it_hint=tmp_ps.g_s_index().end();
        BOOST_FOREACH(term_type t,derived_cast->g_s_index())
        {
          tmp_term=t;
          tmp_term.s_cf()->mult_by_self(cf);
          it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
        }
        derived_cast->swap(tmp_ps);
      }
      template <class DerivedPs2>
        void multiply_terms(const DerivedPs2 &ps2, DerivedPs &retval) const
      {
// TODO: use typedeffed "const_iterator" here instead.
        typedef typename DerivedPs::ancestor::it_s_index it_s_index;
        typedef typename DerivedPs::ancestor::term_type term_type;
        typedef typename DerivedPs2::ancestor::it_s_index it_s_index2;
        typedef typename DerivedPs2::ancestor::term_type term_type2;
        typedef typename DerivedPs::ancestor::cf_type cf_type;
        typedef typename DerivedPs::ancestor::trig_type trig_type;
        typedef typename DerivedPs::ancestor::allocator_type allocator_type;
        typedef light_term<cf_type,trig_type> light_term_type;
        typedef boost::tuple<light_term_type &, light_term_type &> light_term_pair;
        typedef mult_hash<light_term_type,light_term_hasher,
          std::equal_to<light_term_type>,allocator_type,true> m_hash;
        typedef typename m_hash::iterator m_hash_iterator;
        typedef typename m_hash::point_iterator m_hash_point_iterator;
        const DerivedPs *derived_cast=static_cast<DerivedPs const *>(this);
        m_hash hm((derived_cast->length()*ps2.length())/100);
        m_hash_point_iterator hm_p_it;
        const double Delta=derived_cast->g_norm()*ps2.g_norm()*settings_manager::prec(),
          Delta_threshold=Delta/(2*derived_cast->length()*ps2.length());
        size_t n=0;
// NOTE: at this point retval's width() is greater or equal to _both_ this
// and ps2. It's the max width indeed.
        p_assert(math::max(derived_cast->trig_width(),ps2.trig_width())==retval.trig_width());
        light_term_type tmp1, tmp2;
        tmp1.s_trig()->increase_size(retval.trig_width());
        tmp2.s_trig()->increase_size(retval.trig_width());
        light_term_pair term_pair(tmp1,tmp2);
// Cache all pointers to the terms of this and ps2 in vectors.
        std::valarray<term_type const *> v_p1;
        std::valarray<term_type2 const *> v_p2;
        utils::array_pointer(*derived_cast,v_p1);
        utils::array_pointer(ps2,v_p2);
        double norm1;
        size_t i,j;
        const size_t l1=v_p1.size();
        const size_t l2=v_p2.size();
        p_assert(l1 == ps1.length());
        p_assert(l2 == ps2.length());
// Try to build the generalized lexicographic representation.
        series_gl_rep <std::valarray<term_type const *>,
          std::valarray<term_type2 const *>,
          typename DerivedPs::ancestor::trig_type::value_type>
          glr(v_p1,v_p2,derived_cast->trig_width());
        if (glr.is_viable())
        {
          std::cout << "Can do fast" << '\n';
        }
        else
        {
          std::cout << "Can't do fast" << '\n';
// v_it2[0] is legal because we checked for ps2's size.
          const double norm2_i=v_p2[0]->g_cf()->norm(ps2.cf_s_vec());
// Cache some pointers.
          trig_type const *t0=term_pair.template get<0>().g_trig(), *t1=term_pair.template get<1>().g_trig();
          cf_type const *c0=term_pair.template get<0>().g_cf(), *c1=term_pair.template get<1>().g_cf();
          light_term_type *term0=&(term_pair.template get<0>()), *term1=&(term_pair.template get<1>());
          for (i=0;i<l1;++i)
          {
            norm1=v_p1[i]->g_cf()->norm(derived_cast->cf_s_vec());
            if ((norm1*norm2_i)/2<Delta_threshold)
            {
              break;
            }
            for (j=0;j<l2;++j)
            {
// We are going to calculate a term's norm twice... We need to profile
// this at a later stage and see if it is worth to store the norm inside
// the term.
              if ((norm1*v_p2[j]->g_cf()->norm(ps2.cf_s_vec()))/2<Delta_threshold)
              {
                break;
              }
              term_by_term_multiplication(*v_p1[i],*v_p2[j],term_pair);
// Before insertion we change the sign of trigonometric parts if necessary.
// This way we won't do a copy inside the insertion function.
              if (t0->sign()<0)
              {
                term0->invert_trig_args();
              }
              if (t1->sign()<0)
              {
                term1->invert_trig_args();
              }
              hm_p_it=hm.find(*term0);
              if (hm_p_it == hm.end())
              {
                hm.insert(*term0);
              }
              else
              {
                hm_p_it->cf+=*c0;
              }
              hm_p_it=hm.find(*term1);
              if (hm_p_it == hm.end())
              {
                hm.insert(*term1);
              }
              else
              {
                hm_p_it->cf+=*c1;
              }
              ++n;
            }
          }
          const m_hash_iterator hm_it_f=hm.end();
          for (m_hash_iterator hm_it=hm.begin();hm_it!=hm_it_f;++hm_it)
          {
            retval.insert_no_sign_check(term_type(hm_it->cf,hm_it->trig));
          }
//retval.cumulative_crop(Delta);
          std::cout << "w/o trunc=" << derived_cast->length()*ps2.length() << "\tw/ trunc=" << n << std::endl;
          std::cout << "Out length=" << retval.length() << std::endl;
        }
      }
    private:
      // Boilerplate for series multiplication.
      struct light_term_hasher
      {
        size_t operator()(const light_term<typename DerivedPs::cf_type,typename DerivedPs::trig_type> &t) const
        {
          return t.trig.hasher();
        }
      };
      template <class T,class U, class V>
        static void term_by_term_multiplication(const T &t1, const U &t2, V &term_pair)
      {
        typedef typename DerivedPs::ancestor::cf_type cf_type;
        cf_type new_c=*t1.g_cf();
        new_c.mult_by_self(*t2.g_cf());
        new_c/=2;
        DerivedPs::term_by_term_multiplication_trig(t1,t2,term_pair,new_c);
      }
  };
}

#endif
