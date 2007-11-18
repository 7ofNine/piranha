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

#ifndef PIRANHA_TRIGONOMETRIC_TOOLBOX_H
#define PIRANHA_TRIGONOMETRIC_TOOLBOX_H

namespace piranha
{
/// Toolbox for trigonometric operations common to real and complex series.
  template <class DerivedPs>
    class base_trigonometric_toolbox
  {
    public:
      template <class real_DerivedPs>
        void add_ps_to_arg(trig_size_t sym_index, const real_DerivedPs &p)
      {
        typedef std::complex<real_DerivedPs> complex_ps;
        typedef typename DerivedPs::ancestor::it_s_index it_s_index;
        DerivedPs *derived_cast=static_cast<DerivedPs *>(this);
        if (sym_index>=derived_cast->trig_width())
        {
          std::cout << "Invalid index in 'basic_add_ps_to_arg', returning same series." << std::endl;
          return;
        }
        DerivedPs retval;
// Model retval after this.
        action_assert(retval.merge_args(*derived_cast));
// Import also linear arguments.
        retval.lin_args()=derived_cast->lin_args();
// Merge with other series. If we don't succeed just return *this.
        if (!retval.merge_args(p))
        {
          return;
        }
        int16 tmp_mult;
        const it_s_index it_f=derived_cast->end();
        for (it_s_index it=derived_cast->begin();it!=it_f;++it)
        {
          tmp_mult=it->g_trig()->at(sym_index);
// If the symbol's multiplier is zero we simply insert the term.
          if (tmp_mult==0)
          {
            retval.insert(*it);
          }
          else
          {
            complex_ps psc=(p*tmp_mult).complexp();
            real_DerivedPs cosp=psc.real(), sinp=psc.imag();
            DerivedPs tmp1, tmp2;
            action_assert(tmp1.merge_args(*derived_cast));
            tmp1.insert(*it);
            action_assert(tmp2.merge_args(*derived_cast));
            tmp2.insert(*it);
            switch (it->g_flavour())
            {
              case true:
// Change tmp2's flavour.
                tmp2.set_flavour(false);
                retval+=(tmp1*=cosp);
                retval-=(tmp2*=sinp);
                break;
              case false:
                tmp2.set_flavour(true);
                retval+=(tmp1*=cosp);
                retval+=(tmp2*=sinp);
            }
          }
        }
// Take care of lin_args.
        if (derived_cast->lin_args()[sym_index]!=0)
        {
          retval+=(p*derived_cast->
            lin_args()[sym_index]);
        }
        derived_cast->swap(retval);
        return;
      }
      template <class real_DerivedPs>
        void add_ps_to_arg(const std::string &name, const real_DerivedPs &p)
      {
        add_ps_to_arg(static_cast<DerivedPs *>(this)->trig_index(name),p);
      }
  };

/// Toolbox for trigonometric operations.
  template <class DerivedPs>
    class trigonometric_toolbox:public base_trigonometric_toolbox<DerivedPs>
  {
    public:
      typedef std::complex<DerivedPs> complex_ps;
// Maths
      complex_ps complexp() const
      {
        typedef typename complex_ps::ancestor::term_type complex_term_type;
        typedef typename complex_ps::ancestor::cf_type complex_cf_type;
        typedef typename DerivedPs::ancestor::r_it_s_index real_r_it_s_index;
        const DerivedPs *derived_cast=static_cast<DerivedPs const *>(this);
        complex_ps retval;
        action_assert(retval.merge_args(*derived_cast));
        p_assert(retval.trig_width()==derived_cast->trig_width());
        retval.insert(complex_term_type(complex_cf_type(1.)));
        real_r_it_s_index it=derived_cast->g_series_set()->rbegin();
        for (;it!=derived_cast->g_series_set()->rend();++it)
        {
          retval*=jacangdev(it);
          std::cout << retval.length() << std::endl;
        }
        std::cout << "Final size=" << retval.length() << std::endl;
        if (!math::is_zero_vec(derived_cast->lin_args()))
        {
          retval*=complexp_linargs();
        }
        return retval;
      }
/// Calculate cosine of series.
      DerivedPs cosine() const
      {
        return complexp().real();
      }
/// Calculate sine of series.
      DerivedPs sine() const
      {
        return complexp().imag();
      }
    private:
// Maths
/// Complex exponential of a vector of integers.
/**
 * This will generate a single-term Poisson series with unity coefficient.
 */
      complex_ps complexp_linargs() const
      {
        typedef typename complex_ps::ancestor::term_type complex_term_type;
        typedef typename complex_ps::ancestor::cf_type complex_cf_type;
        const DerivedPs *derived_cast=static_cast<DerivedPs const *>(this);
        complex_ps retval;
        action_assert(retval.merge_args(*derived_cast));
        p_assert(retval.trig_width()==derived_cast->trig_width());
        complex_term_type term1(complex_cf_type(1)), term2(complex_cf_type(0,1));
        term2.s_flavour()=false;
        term1.s_trig()->increase_size(retval.trig_width());
        term2.s_trig()->increase_size(retval.trig_width());
        term1.s_trig()->assign_mult_vector(derived_cast->lin_args());
        term2.s_trig()->assign_mult_vector(derived_cast->lin_args());
        retval.insert(term1);
        retval.insert(term2);
        return retval;
      }
      template <class real_cf_type>
        void jaccosRecf(unsigned int i, const real_cf_type &cf, real_cf_type &retval) const
      {
        retval=cf.besselJ(2*i,static_cast<DerivedPs const *>(this)->cf_s_vec());
        retval*=(2-math::Kdelta((unsigned int)0,i))*math::cs_phase(i);
      }
      template <class real_cf_type>
        void jaccosImcf(unsigned int i, const real_cf_type &cf, real_cf_type &retval) const
      {
        retval=cf.besselJ(2*i+1,static_cast<DerivedPs const *>(this)->cf_s_vec());
        retval*=2*math::cs_phase(i);
      }
      template <class real_cf_type>
        void jacsinRecf(unsigned int i, const real_cf_type &cf, real_cf_type &retval) const
      {
        retval=cf.besselJ(2*i,static_cast<DerivedPs const *>(this)->cf_s_vec());
        retval*=(2-math::Kdelta((unsigned int)0,i));
      }
      template <class real_cf_type>
        void jacsinImcf(unsigned int i, const real_cf_type &cf, real_cf_type &retval) const
      {
        retval=cf.besselJ(2*i+1,static_cast<DerivedPs const *>(this)->cf_s_vec());
        retval*=2;
      }
/// Jacobi-Anger development of a term.
/**
 * Mathematically speaking this development converges gracely for series whose coefficients are
 * <=1. It should not be a problem with ephemerides of celestial bodies - the 1st order angular
 * perturbation(s) must be less than 1rad/timeunit - maybe then it is enough to change timescale?
 * Or maybe we can split up the biggest terms, do separate complexps and multiply at the end
 * and assemble the result....
 */
// This has to be templated this way because we don't know, during multiple inheritance,
// about typedefs in DerivedPs class. Fortunately the problem is not serious here, since this
// function is private and when it is called the compiler determines automatically the iterator
// type from the context.
      template <class Iterator>
        complex_ps jacangdev(Iterator it) const
      {
        typedef typename complex_ps::ancestor::term_type complex_term_type;
        typedef typename complex_ps::ancestor::cf_type complex_cf_type;
        typedef typename DerivedPs::ancestor::cf_type real_cf_type;
        unsigned int i;
        complex_ps retval;
        action_assert(retval.merge_args(*static_cast<DerivedPs const *>(this)));
        p_assert(retval.trig_width()==static_cast<DerivedPs const *>(this)->trig_width());
        real_cf_type _cf=*it->g_cf(), tmp;
        complex_term_type term1, term2;
        if (it->g_flavour())
        {
          for (i=0;i<settings_manager::jacang_lim();++i)
          {
            jaccosRecf<real_cf_type>(i,_cf,tmp);
            term1.s_cf()->set_real(tmp);
            *term1.s_trig()=*it->g_trig();
            *term1.s_trig()*=(i<<1);
            term1.s_flavour()=true;
            retval.insert(term1);
            jaccosImcf<real_cf_type>(i,_cf,tmp);
            term2.s_cf()->set_imag(tmp);
            *term2.s_trig()=*it->g_trig();
            *term2.s_trig()*=((i<<1)+1);
            term1.s_flavour()=true;
            retval.insert(term2);
          }
        }
        else
        {
          for (i=0;i<settings_manager::jacang_lim();++i)
          {
            jacsinRecf(i,_cf,tmp);
            term1.s_cf()->set_real(tmp);
            *term1.s_trig()=*it->g_trig();
            *term1.s_trig()*=(i<<1);
            term1.s_flavour()=true;
            retval.insert(term1);
            jacsinImcf(i,_cf,tmp);
            term2.s_cf()->set_imag(tmp);
            *term2.s_trig()=*it->g_trig();
            *term2.s_trig()*=((i<<1)+1);
            term2.s_flavour()=false;
            retval.insert(term2);
          }
        }
        return retval;
      }
  };

/// Toolbox for trigonometric operations, specialization for complex series.
  template <>
    template <class DerivedPs>
    class trigonometric_toolbox<std::complex<DerivedPs> >:public base_trigonometric_toolbox<std::complex<DerivedPs> >
  {};
}
#endif
