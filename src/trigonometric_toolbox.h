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
  template <class Derived, class real_Derived>
    class common_trig_toolbox
  {
    public:
      typedef std::complex<real_Derived> complex_ps;
      void add_ps_to_arg(trig_size_t sym_index, const real_Derived &p)
      {
        typedef typename Derived::ancestor::it_s_index it_s_index;
        if (sym_index>=static_cast<Derived *>(this)->
          trig_width())
        {
          std::cout << "Invalid index in 'basic_add_ps_to_arg', returning same series." << std::endl;
          return;
        }
        Derived retval;
// Model retval after this.
        p_assert(retval.merge_args(*static_cast<Derived *>(this)));
// Import also linear arguments.
        retval.lin_args()=static_cast<Derived *>(this)->lin_args();
// Merge with other series. If we don't succeed just return *this.
        if (!retval.merge_args(p))
        {
          return;
        }
        mult_t tmp_mult;
// TODO: cache end() here.
        for (it_s_index
          it=static_cast<Derived *>(this)->begin();it!=static_cast<Derived *>(this)->end();++it)
        {
          tmp_mult=it->trig_args().multiplier(sym_index);
// If the symbol's multiplier is zero we simply insert the term.
          if (tmp_mult==0)
          {
            retval.insert(*it);
          }
          else
          {
            complex_ps psc=(p*tmp_mult).complexp();
            real_Derived cosp=psc.real(), sinp=psc.imag();
            Derived tmp1, tmp2;
            p_assert(tmp1.merge_args(*static_cast<Derived *>(this)));
            tmp1.insert(*it);
            p_assert(tmp2.merge_args(*static_cast<Derived *>(this)));
            tmp2.insert(*it);
            switch (it->flavour())
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
        if (static_cast<Derived *>(this)->lin_args()[sym_index]!=0)
        {
          retval+=(p*static_cast<Derived *>(this)->
            lin_args()[sym_index]);
        }
        static_cast<Derived *>(this)->swap(retval);
        return;
      }
      void add_ps_to_arg(const std::string &name, const real_Derived &p)
      {
        add_ps_to_arg(static_cast<Derived *>(this)->trig_index(name),p);
      }
  };

/// Toolbox for trigonometric operations specific to real series.
  template <class Derived>
    class real_trig_toolbox
  {
    public:
      typedef std::complex<Derived> complex_ps;
// Maths
      complex_ps complexp() const
      {
        typedef typename complex_ps::ancestor::term_type complex_term_type;
        typedef typename complex_ps::ancestor::cf_type complex_cf_type;
        typedef typename Derived::ancestor::r_it_s_index real_r_it_s_index;
        complex_ps retval;
        p_assert(retval.merge_args(*static_cast<Derived const *>(this)));
        p_assert(retval.trig_width()==static_cast<Derived const *>(this)->trig_width());
        retval.insert(complex_term_type(complex_cf_type(1.)));
        real_r_it_s_index it=static_cast<Derived const *>(this)->set().rbegin();
        for (;it!=static_cast<Derived const *>(this)->set().rend();++it)
        {
          retval*=jacangdev(it);
          std::cout << retval.length() << std::endl;
        }
        std::cout << "Final size=" << retval.length() << std::endl;
        if (!is_zero_vec(static_cast<Derived const *>(this)->lin_args()))
        {
          retval*=complexp_linargs()
            ;
        }
        return retval;
      }
/// Complex multiangle.
      complex_ps complex_multiangle(trig_size_t pos, mult_t n) const
      {
        typedef typename complex_ps::ancestor::term_type complex_term_type;
        typedef typename complex_ps::ancestor::cf_type complex_cf_type;
        if (pos>=static_cast<Derived const *>(this)->
          trig_width())
        {
          std::cout << "Error: requested multiangle complex with invalid integer position" <<
            std::endl;
          std::exit(1);
        }
        complex_ps retval;
        retval.merge_args(*static_cast<Derived const *>(this));
        complex_term_type in_term(complex_cf_type(1),true);
        in_term.trig_args().increase_size(static_cast<Derived const *>(this)->trig_width());
        in_term.trig_args().insert(pos,n);
        retval.insert(in_term);
        in_term.flavour()=false;
        retval.insert(in_term);
        return retval;
      }
/// Calculate cosine of series.
      Derived cosine() const
      {
        return complexp().real();
      }
/// Calculate sine of series.
      Derived sine() const
      {
        return complexp().imag();
      }
//|-------------------------|
//|Protected Specializations|
//|-------------------------|
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
        complex_ps retval;
        p_assert(retval.merge_args(*static_cast<Derived const *>(this)));
        p_assert(retval.trig_width()==static_cast<Derived const *>(this)->trig_width());
        complex_term_type term1(complex_cf_type(1),true),
          term2(complex_cf_type(0,1),false);
        term1.trig_args().increase_size(retval.trig_width());
        term2.trig_args().increase_size(retval.trig_width());
        for (unsigned int j=0;
          j<static_cast<Derived const *>
          (this)->lin_args().size();
          ++j)
        {
          term1.trig_args().insert(j,static_cast<Derived const *>(this)->lin_args()[j]);
          term2.trig_args().insert(j,static_cast<Derived const *>(this)->lin_args()[j]);
        }
        retval.insert(term1);
        retval.insert(term2);
        return retval;
      }
/// Jacobi-Anger development of a term.
/**
 * Mathematically speaking this development converges gracely for series whose coefficients are
 * <=1. It should not be a problem with ephemerides of celestial bodies - the 1st order angular
 * perturbation(s) must be less than 1rad/timeunit - maybe then it is enough to change timescale?
 * Or maybe we can split up the biggest terms, do separate complexps and multiply at the end
 * and assemble the result....
 * If we start using these besselJ inside symbolic manipulation we need to turn the macro below
 * into proper functions.
 */
#define __jaccosRecf(i,cf) (cf.besselJ(2*i,static_cast<Derived const *>(this)->cf_s_vec())*(2-math::Kdelta((unsigned int)0,i))*math::cs_phase(i))
#define __jaccosImcf(i,cf) (cf.besselJ(2*i+1,static_cast<Derived const *>(this)->cf_s_vec())*2*math::cs_phase(i))

#define __jacsinRecf(i,cf) (cf.besselJ(2*i,static_cast<Derived const *>(this)->cf_s_vec())*(2-math::Kdelta((unsigned int)0,i)))
#define __jacsinImcf(i,cf) (cf.besselJ(2*i+1,static_cast<Derived const *>(this)->cf_s_vec())*2)
// This has to be templated this way because we don't know, during multiple inheritance,
// about typedefs in Derived class. Fortunately the problem is not serious here, since this
// function is private and when it is called the compiler determines automatically the iterator
// type from the context.
      template <class Iterator>
        complex_ps jacangdev(Iterator it) const
      {
        typedef typename complex_ps::ancestor::term_type complex_term_type;
        typedef typename complex_ps::ancestor::cf_type complex_cf_type;
        typedef typename Derived::ancestor::cf_type real_cf_type;
        unsigned int i, w=static_cast<Derived const *>(this)->trig_width();
        complex_ps retval;
        p_assert(retval.merge_args(*static_cast<Derived const *>(this)));
        p_assert(retval.trig_width()==static_cast<Derived const *>(this)->trig_width());
        real_cf_type _cf=it->c();
        if (it->flavour()
          )
        {
          complex_term_type term1, term2;
          term1.trig_args().increase_size(w);
          term2.trig_args().increase_size(w);
          for (i=0;i<settings_manager::jacang_lim();++i)
          {
            term1.c().set_real(__jaccosRecf(i,_cf));
            term1.trig_args()=it->trig_args();
            term1.trig_args()*=(i<<1);
            retval.insert(term1);
            term2.c().set_imag(__jaccosImcf(i,_cf));
            term2.trig_args()=it->trig_args();
            term2.trig_args()*=((i<<1)+1);
            retval.insert(term2);
          }
        }
        else
        {
          complex_term_type term1, term2(complex_cf_type(0),false);
          term1.trig_args().increase_size(w);
          term2.trig_args().increase_size(w);
          for (i=0;i<settings_manager::jacang_lim();++i)
          {
            term1.c().set_real(__jacsinRecf(i,_cf));
            term1.trig_args()=it->trig_args();
            term1.trig_args()*=(i<<1);
            retval.insert(term1);
            term2.c().set_imag(__jacsinImcf(i,_cf));
            term2.trig_args()=it->trig_args();
            term2.trig_args()*=((i<<1)+1);
            retval.insert(term2);
          }
        }
        return retval;
      }
#undef __jaccosRecf
#undef __jaccosImcf
#undef __jacsinRecf
#undef __jacsinImcf

  };
}
#endif
