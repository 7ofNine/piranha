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

#ifndef PIRANHA_NUMERICAL_ELEMENTARY_MATH_TOOLBOX_H
#define PIRANHA_NUMERICAL_ELEMENTARY_MATH_TOOLBOX_H

namespace piranha
{
/// Elementary math toolbox for numerical series.
/**
 * Implements series multiplication and truncation according to norms.
 */
  template <class Derived>
    class numerical_elementary_math_toolbox
  {
    public:
      template <class Derived2>
        void series_multiplication(const Derived2 &ps2)
      {
        Derived *derived_cast=static_cast<Derived *>(this);
        Derived retval;
// NOTE: optimize in the case one series is a c value?
// If one length is zero do not do anything
        if (derived_cast->length()!=0 && ps2.length()!=0)
        {
          if (!is_zero_vec(derived_cast->lin_args())||!is_zero_vec(ps2.lin_args()))
          {
            std::cout << "Non-zero linargs!" << std::endl;
            std::exit(1);
          }
          if (!derived_cast->merge_args(ps2))
          {
            std::cout << "args are not compatible, returning self." << std::endl;
            std::exit(1);
            return;
          }
          p_assert(retval.merge_args(*derived_cast));
          const double Delta=derived_cast->norm()*ps2.norm()*settings_manager::prec();
          multiply_terms(ps2,retval,Delta);
        }
        else
        {
          std::cout << "Zero stuff" << std::endl;
        }
        derived_cast->swap(retval);
      }
    private:
      template <class Derived2>
        void multiply_terms(const Derived2 &ps2, Derived &retval, const double &Delta) const
      {
        typedef typename Derived::ancestor::it_s_index it_s_index;
        typedef typename Derived2::ancestor::it_s_index it_s_index2;
        typedef typename Derived::ancestor::term_type term_type;
        const Derived *derived_cast=static_cast<Derived const *>(this);
        const double Delta_threshold=Delta/(2*derived_cast->length()*ps2.length());
        size_t n=0;
// NOTE: at this point retval's width() is greater or equal to _both_ this
// and ps2. It's the max width indeed.
        p_assert(math::max(derived_cast->trig_width(),ps2.trig_width())==retval.trig_width());
        term_type tmp1, tmp2;
        boost::tuple<term_type &, term_type &> term_pair(tmp1,tmp2);
        const it_s_index it1_f=derived_cast->s_index().end();
        const it_s_index2 it2_f=ps2.s_index().end();
        it_s_index2 it2;
        it_s_index it1, it_hint=retval.s_index().end();
        double norm1;
        for (it1=derived_cast->s_index().begin();it1!=it1_f;++it1)
        {
          it2=ps2.s_index().begin();
          norm1=it1->norm(derived_cast->cf_s_vec());
          if ((norm1*it2->norm(ps2.cf_s_vec()))/2<Delta_threshold)
          {
            break;
          }
          for (;it2!=it2_f;++it2)
          {
// We are going to calculate a term's norm twice... We need to profile
// this at a later stage and see if it is worth to store the norm inside
// the term.
            if ((norm1*it2->norm(ps2.cf_s_vec()))/2<Delta_threshold)
            {
              break;
            }
            term_by_term_multiplication(*it1,*it2,term_pair);
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
            ++n;
          }
        }
//retval.cumulative_crop(Delta);
        std::cout << "w/o trunc=" << derived_cast->length()*ps2.length() << "\tw/ trunc=" << n << std::endl;
        std::cout << "Out length=" << retval.length() << std::endl;
      }
    template <class T,class U>
      void term_by_term_multiplication(const T &t1, const U &t2, boost::tuple<T &,T &> &term_pair) const
    {
      typedef typename Derived::ancestor::cf_type cf_type;
      cf_type new_c=t1.g_cf();
      new_c*=t2.g_cf();
      new_c/=2;
      static_cast<Derived const *>(this)->term_by_term_multiplication_trig(t1,t2,term_pair,new_c);
    }
  };
}
#endif
