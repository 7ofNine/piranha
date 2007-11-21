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

#ifndef FOURIER_MULTIPLICATION_TOOLBOX_H
#define FOURIER_MULTIPLICATION_TOOLBOX_H

#include <boost/foreach.hpp>
#include <cstring> // For memset.

#include "../bits/buffer.h"
#include "../bits/config.h" // For selection of temporary hash container for multiplication.
#include "../bits/light_term.h"
#include "../bits/norm_truncation.h"
#include "../bits/pseries_gl_rep.h"

namespace piranha
{
/// Multiplication toolbox for Fourier series.
/**
 * Series multiplication used a norm-based truncation strategy.
 */
  template <class DerivedPs>
    class fourier_multiplication_toolbox:public norm_truncation
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
// DerivedPs typedefs.
        typedef typename DerivedPs::ancestor::it_s_index it_s_index;
        typedef typename DerivedPs::ancestor::term_type term_type;
        typedef typename DerivedPs2::ancestor::it_s_index it_s_index2;
        typedef typename DerivedPs2::ancestor::term_type term_type2;
        typedef typename DerivedPs::ancestor::cf_type cf_type;
        typedef typename DerivedPs::ancestor::trig_type trig_type;
        typedef typename DerivedPs::ancestor::allocator_type allocator_type;
// Light_term typedefs.
        typedef light_term<cf_type,trig_type> light_term_type;
        typedef boost::tuple<light_term_type &, light_term_type &> light_term_pair;
// Mult_hash typedefs.
        typedef mult_hash<light_term_type,typename light_term_type::hasher,
          std::equal_to<light_term_type>,allocator_type,true> m_hash;
        typedef typename m_hash::iterator m_hash_iterator;
        typedef typename m_hash::point_iterator m_hash_point_iterator;
// Pseries_gl_rep typedefs.
        typedef pseries_gl_rep<DerivedPs,DerivedPs2> glr_type;
        const DerivedPs *derived_cast=static_cast<DerivedPs const *>(this);
        const size_t l1=derived_cast->length(), l2=ps2.length();
        const double Delta=derived_cast->g_norm()*ps2.g_norm()*get_truncation(),
          Delta_threshold=Delta/(2*l1*l2);
        p_assert(math::max(derived_cast->trig_width(),ps2.trig_width()) == retval.trig_width());
        double norm1;
        size_t i,j;
// ps2.begin() is legal because we checked for ps2's size.
        const double norm2_i=ps2.begin()->g_cf()->norm(ps2.cf_s_vec());
// Build the generalized lexicographic representation.
        glr_type glr(*derived_cast,ps2);
        if (glr.is_viable())
        {
          const max_fast_int h_min = glr.g_h_min(), h_max = glr.g_h_max();
          const max_fast_int h_card = (h_max - h_min +1);
          p_assert(h_card >= 0);
          std::cout << "h_card: " << h_card << '\n';
          std::cout << "h_minmax: " << h_min << ',' << h_max << '\n';
          typedef std::pair<cf_type,bool> cf_bool;
          const double load_factor = ((double)l1*l2)/h_card;
          std::cout << "Load factor: " << load_factor << '\n';
// TODO: hard-wire this for now, we have to study it a bit.
#define _MAX_LOAD_FACTOR (1.)
          if (load_factor < _MAX_LOAD_FACTOR)
          {
            std::cout << "Load factor is too small, will avoid coded vector arithmetics." << std::endl;
          }
          if ((size_t)(h_card<<1) <= buffer::n_elements<cf_bool>() and load_factor >= _MAX_LOAD_FACTOR)
#undef _MAX_LOAD_FACTOR
          {
            std::cout << "Can do vector coded arithmetics." << '\n';
            coded_vector_mult(h_card,h_min,h_max,norm2_i,Delta_threshold,glr,retval,l1,l2,ps2);
          }
          else
          {
            std::cout << "Can do hash coded arithmetics." << '\n';
            coded_hash_mult(glr,l1,l2,retval,norm2_i,Delta_threshold,ps2);
          }
        }
        else
        {
          std::cout << "Can't do fast" << '\n';
          size_t n=0;
// NOTE: at this point retval's width() is greater or equal to _both_ this
// and ps2. It's the max width indeed.
          light_term_type tmp1, tmp2;
          tmp1.s_trig()->increase_size(retval.trig_width());
          tmp2.s_trig()->increase_size(retval.trig_width());
          light_term_pair term_pair(tmp1,tmp2);
// Cache all pointers to the terms of this and ps2 in vectors.
          std::valarray<term_type const *> v_p1;
          std::valarray<term_type2 const *> v_p2;
          utils::array_pointer(*derived_cast,v_p1);
          utils::array_pointer(ps2,v_p2);
          p_assert(v_p1.size() == derived_cast->length());
          p_assert(v_p2.size() == ps2.length());
// Prepare the structure for multiplication.
          m_hash hm((derived_cast->length()*ps2.length())/100);
          m_hash_point_iterator hm_p_it;
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
// TODO possible optimization: introduce destructive term ctor (e.g., swap array content instead of copying it)?
            retval.insert_no_sign_check(term_type(hm_it->cf,hm_it->trig));
          }
//retval.cumulative_crop(Delta);
          std::cout << "w/o trunc=" << derived_cast->length()*ps2.length() << "\tw/ trunc=" << n << std::endl;
          std::cout << "Out length=" << retval.length() << std::endl;
        }
      }
    private:
      template <class T,class U, class V>
        static void term_by_term_multiplication(const T &t1, const U &t2, V &term_pair)
      {
        typedef typename DerivedPs::ancestor::cf_type cf_type;
        cf_type new_c=*t1.g_cf();
        new_c.mult_by_self(*t2.g_cf());
        new_c/=2;
        DerivedPs::term_by_term_multiplication_trig(t1,t2,term_pair,new_c);
      }
      template <class Glr, class Retval, class DerivedPs2>
        void coded_vector_mult(const max_fast_int &h_card, const max_fast_int &h_min, const max_fast_int &h_max,
        const double &norm2_i, const double &Delta_threshold, const Glr &glr,
        Retval &retval, const size_t &l1, const size_t &l2, const DerivedPs2 &ps2) const
      {
        typedef typename DerivedPs::ancestor::cf_type cf_type;
        typedef typename DerivedPs::ancestor::trig_type::value_type mult_type;
        typedef typename DerivedPs::ancestor::term_type term_type;
        typedef Glr glr_type;
        typedef typename glr_type::coded_series_type1 cs_type1;
        typedef typename glr_type::coded_series_type2 cs_type2;
        typedef std::pair<cf_type,bool> cf_bool;
        const DerivedPs *derived_cast=static_cast<DerivedPs const *>(this);
        const cs_type1 &cs1 = glr.g1();
        const cs_type2 &cs2 = glr.g2();
        cf_bool *code_vector_cos = buffer::head<cf_bool>(),
        *code_vector_sin = code_vector_cos + h_card;
// Reset memory area.
        memset(code_vector_cos,0,sizeof(cf_bool)*h_card);
        memset(code_vector_sin,0,sizeof(cf_bool)*h_card);
        cf_bool *s_point_cos = code_vector_cos - h_min,
          *s_point_sin = code_vector_sin - h_min;
        cf_type tmp_cf;
        size_t i, j;
        double norm1;
        for (i=0;i<l1;++i)
        {
          norm1=cs1[i].cf.norm(derived_cast->cf_s_vec());
          if ((norm1*norm2_i)/2<Delta_threshold)
          {
            break;
          }
          for (j=0;j<l2;++j)
          {
            if ((norm1*cs2[j].cf.norm(ps2.cf_s_vec()))/2<Delta_threshold)
            {
              break;
            }
            tmp_cf = cs1[i].cf;
            tmp_cf.mult_by_self(cs2[j].cf);
            tmp_cf/=2;
            const max_fast_int tmp_index_plus = cs1[i].code + cs2[j].code,
              tmp_index_minus = cs1[i].code - cs2[j].code;
            if (cs1[i].flavour == cs2[j].flavour)
            {
              if (cs1[i].flavour)
              {
                s_point_cos[tmp_index_minus].first += tmp_cf;
                s_point_cos[tmp_index_minus].second = true;
                s_point_cos[tmp_index_plus].first += tmp_cf;
                s_point_cos[tmp_index_plus].second = true;
              }
              else
              {
                s_point_cos[tmp_index_minus].first += tmp_cf;
                s_point_cos[tmp_index_minus].second = true;
                s_point_cos[tmp_index_plus].first -= tmp_cf;
                s_point_cos[tmp_index_plus].second = true;
              }
            }
            else
            {
              if (cs1[i].flavour)
              {
                s_point_sin[tmp_index_minus].first -= tmp_cf;
                s_point_sin[tmp_index_minus].second = true;
                s_point_sin[tmp_index_plus].first += tmp_cf;
                s_point_sin[tmp_index_plus].second = true;
              }
              else
              {
                s_point_sin[tmp_index_minus].first += tmp_cf;
                s_point_sin[tmp_index_minus].second = true;
                s_point_sin[tmp_index_plus].first += tmp_cf;
                s_point_sin[tmp_index_plus].second = true;
              }
            }
          }
        }
        term_type tmp_term;
        tmp_term.s_trig()->increase_size(derived_cast->trig_width());
        std::valarray<mult_type> tmp_array(derived_cast->trig_width());
        max_fast_int k;
        for (k=h_min;k<=h_max;++k)
        {
          if (unlikely(s_point_cos[k].second))
          {
            *tmp_term.s_cf() = s_point_cos[k].first;
            glr.decode_multiindex(k,tmp_array);
            tmp_term.s_trig()->assign_mult_vector(tmp_array);
            tmp_term.s_trig()->s_flavour()=true;
            retval.insert(tmp_term);
          }
        }
        for (k=h_min;k<=h_max;++k)
        {
          if (unlikely(s_point_sin[k].second))
          {
            *tmp_term.s_cf() = s_point_sin[k].first;
            glr.decode_multiindex(k,tmp_array);
            tmp_term.s_trig()->assign_mult_vector(tmp_array);
            tmp_term.s_trig()->s_flavour()=false;
            retval.insert(tmp_term);
          }
        }
        std::cout << "Out length=" << retval.length() << std::endl;
      }
      template <class Glr, class Retval, class DerivedPs2>
        void coded_hash_mult(const Glr &glr, const size_t &l1, const size_t &l2,
        Retval &retval, const double &norm2_i, const double &Delta_threshold,
        const DerivedPs2 &ps2) const
      {
        typedef typename DerivedPs::ancestor::term_type term_type;
        typedef typename DerivedPs::ancestor::trig_type::value_type mult_type;
        typedef typename DerivedPs::ancestor::allocator_type allocator_type;
        typedef Glr glr_type;
        typedef typename glr_type::coded_series_type1 cs_type1;
        typedef typename glr_type::coded_series_type2 cs_type2;
        typedef typename glr_type::cct_type1 cct_type;
// Coded mult_hash typedefs.
        typedef mult_hash<cct_type,typename glr_type::cct_hasher,
          typename glr_type::cct_equal_to,
          allocator_type,true> ccm_hash;
        typedef typename ccm_hash::iterator ccm_hash_iterator;
        typedef typename ccm_hash::point_iterator ccm_hash_point_iterator;
        const DerivedPs *derived_cast=static_cast<DerivedPs const *>(this);
        const cs_type1 &cs1 = glr.g1();
        const cs_type2 &cs2 = glr.g2();
        cct_type tmp_term1, tmp_term2;
        ccm_hash cchm_cos((l1*l2)/100),
          cchm_sin((l1*l2)/100);
        ccm_hash_point_iterator cchm_p_it;
        double norm1;
        size_t i, j;
        for (i=0;i<l1;++i)
        {
          norm1=cs1[i].cf.norm(derived_cast->cf_s_vec());
          if ((norm1*norm2_i)/2<Delta_threshold)
          {
            break;
          }
          for (j=0;j<l2;++j)
          {
            if ((norm1*cs2[j].cf.norm(ps2.cf_s_vec()))/2<Delta_threshold)
            {
              break;
            }
            tmp_term1.cf = cs1[i].cf;
            tmp_term2.cf = cs1[i].cf;
            tmp_term1.code = cs1[i].code;
            tmp_term2.code = cs1[i].code;
// For now calculate a+b and a-b.
            tmp_term1.code-=cs2[j].code;
            tmp_term2.code+=cs2[j].code;
// Now the coefficients, all with positive signs for now.
            tmp_term1.cf.mult_by_self(cs2[j].cf);
            tmp_term1.cf/=2;
            tmp_term2.cf=tmp_term1.cf;
// Now fix flavours and coefficient signs.
            if (cs1[i].flavour == cs2[j].flavour)
            {
              if (!(cs1[i].flavour))
              {
                tmp_term2.cf *= -1;
              }
// Insert into cosine container.
              cchm_p_it = cchm_cos.find(tmp_term1);
              if (cchm_p_it == cchm_cos.end())
              {
                cchm_cos.insert(tmp_term1);
              }
              else
              {
                cchm_p_it->cf+=tmp_term1.cf;
              }
              cchm_p_it = cchm_cos.find(tmp_term2);
              if (cchm_p_it == cchm_cos.end())
              {
                cchm_cos.insert(tmp_term2);
              }
              else
              {
                cchm_p_it->cf+=tmp_term2.cf;
              }
            }
            else
            {
              if (cs1[i].flavour)
              {
                tmp_term1.cf *= -1;
              }
// Insert into sine container.
              cchm_p_it = cchm_sin.find(tmp_term1);
              if (cchm_p_it == cchm_sin.end())
              {
                cchm_sin.insert(tmp_term1);
              }
              else
              {
                cchm_p_it->cf+=tmp_term1.cf;
              }
              cchm_p_it = cchm_sin.find(tmp_term2);
              if (cchm_p_it == cchm_sin.end())
              {
                cchm_sin.insert(tmp_term2);
              }
              else
              {
                cchm_p_it->cf+=tmp_term2.cf;
              }
            }
          }
        }
        term_type tmp_term;
        tmp_term.s_trig()->increase_size(derived_cast->trig_width());
        std::valarray<mult_type> tmp_array(derived_cast->trig_width());
        ccm_hash_iterator cchm_it;
        {
          const ccm_hash_iterator cchm_it_f=cchm_cos.end();
          for (cchm_it=cchm_cos.begin();cchm_it!=cchm_it_f;++cchm_it)
          {
            *tmp_term.s_cf() = cchm_it->cf;
            glr.decode_multiindex(cchm_it->code,tmp_array);
            tmp_term.s_trig()->assign_mult_vector(tmp_array);
            tmp_term.s_trig()->s_flavour()=true;
            retval.insert(tmp_term);
          }
        }
        {
          const ccm_hash_iterator cchm_it_f=cchm_sin.end();
          for (cchm_it=cchm_sin.begin();cchm_it!=cchm_it_f;++cchm_it)
          {
            *tmp_term.s_cf() = cchm_it->cf;
            glr.decode_multiindex(cchm_it->code,tmp_array);
            tmp_term.s_trig()->assign_mult_vector(tmp_array);
            tmp_term.s_trig()->s_flavour()=false;
            retval.insert(tmp_term);
          }
        }
        std::cout << "Out length=" << retval.length() << std::endl;
      }
  };
}

#endif
