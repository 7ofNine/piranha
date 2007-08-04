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

#ifndef PIRANHA_SYMBOL_LIMITING_ELEMENTARY_MATH_TOOLBOX_H
#define PIRANHA_SYMBOL_LIMITING_ELEMENTARY_MATH_TOOLBOX_H

namespace piranha
{
/// Elementary math toolbox for symbolic series.
/**
 * Implements series truncation according to symbol limits set in piranha::psymbol_limiter.
 */
  template <class Derived>
    class symbol_limiting_elementary_math_toolbox
  {
    public:
      template <class Derived2>
        void series_multiplication(const Derived2 &ps2)
      {
        typedef typename Derived::ancestor::cf_type cf_type;
        Derived *derived_cast=static_cast<Derived *>(this);
        Derived retval;
        if (derived_cast->series_multiplication_preliminaries(ps2,retval))
        {
          if (derived_cast->series_multiplication_optimize_for_cf(ps2))
          {
            return;
          }
          multiply_terms(ps2,retval);
        }
        derived_cast->swap(retval);
      }
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
        it_s_index it_hint=tmp_ps.s_index().end();
        const it_s_index it_f=derived_cast->s_index().end();
        vec_expo_index_limit limits;
        symbol_limiter::get_limits_index(derived_cast->cf_s_vec(),limits);
        for (it_s_index it=derived_cast->s_index().begin();it!=it_f;++it)
        {
          tmp_term=*it;
          tmp_term.s_cf().mult_by_self(cf,limits);
          it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
        }
        derived_cast->swap(tmp_ps);
      }
    private:
      template <class Derived2>
        void multiply_terms(const Derived2 &ps2, Derived &retval) const
      {
        typedef typename Derived::ancestor::it_s_index it_s_index;
        typedef typename Derived2::ancestor::it_s_index it_s_index2;
        typedef typename Derived::ancestor::term_type term_type;
        const Derived *derived_cast=static_cast<Derived const *>(this);
// NOTE: at this point retval's width() is greater or equal to _both_ this
// and ps2. It's the max width indeed.
        p_assert(math::max(derived_cast->trig_width(),ps2.trig_width())==retval.trig_width());
        term_type tmp1, tmp2;
        boost::tuple<term_type &, term_type &> term_pair(tmp1,tmp2);
        const it_s_index it1_f=derived_cast->s_index().end();
        const it_s_index2 it2_f=ps2.s_index().end();
        it_s_index2 it2;
        const it_s_index it2_i=ps2.s_index().begin();
        it_s_index it1, it_hint=retval.s_index().end();
        vec_expo_index_limit limits;
        symbol_limiter::get_limits_index(derived_cast->cf_s_vec(),limits);
        for (it1=derived_cast->s_index().begin();it1!=it1_f;++it1)
        {
          for (it2=it2_i;it2!=it2_f;++it2)
          {
            term_by_term_multiplication(*it1,*it2,term_pair,limits);
// Before insertion we change the sign of trigonometric parts if necessary.
// This way we won't do a copy inside insertion function.
            if (term_pair.template get
              <0>().g_trig().sign()<0)
            {
              term_pair.template get<0>().invert_trig_args();
            }
            if (term_pair.template get
              <1>().g_trig().sign()<0)
            {
              term_pair.template get<1>().invert_trig_args();
            }
            it_hint=retval.insert(term_pair.template get
              <0>(),true,&it_hint);
            it_hint=retval.insert(term_pair.template get
              <1>(),true,&it_hint);
          }
        }
      }
      template <class T,class U>
        static void term_by_term_multiplication(const T &t1, const U &t2, boost::tuple<T &,T &> &term_pair,
        const vec_expo_index_limit &limits)
      {
        typedef typename Derived::ancestor::cf_type cf_type;
        cf_type new_c=t1.g_cf();
        new_c.mult_by_self(t2.g_cf(),limits);
        new_c/=2;
        Derived::term_by_term_multiplication_trig(t1,t2,term_pair,new_c);
      }
  };
}
#endif
