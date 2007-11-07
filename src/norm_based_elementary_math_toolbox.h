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

#include <boost/foreach.hpp>

#include "config.h"         // For selection of temporary multiplication hash container
#include "light_term.h"

namespace piranha
{
/// Elementary math toolbox for numerical series.
/**
 * Implements series multiplication and truncation according to norms.
 */
  template <class Derived>
    class norm_based_elementary_math_toolbox
  {
    template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class DerivedPs, class Allocator>
      friend class base_pseries;
    private:
      template <class Cf>
        void cf_multiplication(const Cf &cf)
      {
        typedef typename Derived::ancestor::term_type term_type;
        typedef typename Derived::ancestor::it_s_index it_s_index;
        Derived *derived_cast=static_cast<Derived *>(this);
        if (derived_cast->empty())
        {
          return;
        }
        Derived tmp_ps;
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
// Boilerplate for series multiplication.
      struct light_term_hasher
      {
        size_t operator()(const light_term<typename Derived::cf_type,typename Derived::trig_type> &t) const
        {
          return t.trig.hasher();
        }
      };
      template <class Derived2>
        void multiply_terms(const Derived2 &ps2, Derived &retval) const
      {
// TODO: use typedeffed "const_iterator" here instead.
        typedef typename Derived::ancestor::it_s_index it_s_index;
        typedef typename Derived::ancestor::term_type term_type;
        typedef typename Derived2::ancestor::it_s_index it_s_index2;
        typedef typename Derived2::ancestor::term_type term_type2;
        typedef typename Derived::ancestor::cf_type cf_type;
        typedef typename Derived::ancestor::trig_type trig_type;
        typedef typename Derived::ancestor::allocator_type allocator_type;
        typedef light_term<cf_type,trig_type> light_term_type;
        typedef boost::tuple<light_term_type &, light_term_type &> light_term_pair;
        typedef mult_hash<light_term_type,light_term_hasher,
          std::equal_to<light_term_type>,allocator_type,true> m_hash;
        typedef typename m_hash::iterator m_hash_iterator;
        typedef typename m_hash::point_iterator m_hash_point_iterator;
        const Derived *derived_cast=static_cast<Derived const *>(this);
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
// This way we won't do a copy inside insertion function.
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
      template <class T,class U, class V>
        static void term_by_term_multiplication(const T &t1, const U &t2, V &term_pair)
      {
        typedef typename Derived::ancestor::cf_type cf_type;
        cf_type new_c=*t1.g_cf();
        new_c.mult_by_self(*t2.g_cf());
        new_c/=2;
        Derived::term_by_term_multiplication_trig(t1,t2,term_pair,new_c);
      }
  };
}
#endif
