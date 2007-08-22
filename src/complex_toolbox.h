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

#ifndef COMPLEX_TOOLBOX_H
#define COMPLEX_TOOLBOX_H

namespace piranha
{
  template <class real_Derived>
    class complex_toolbox
  {
    protected:
// NOTICE: typedefs regarding Derived type cannot be placed here because when the compiler
// parses this part it does not know enough about Derived yet. They are ok in the body of
// the methods tough, because instantiation happens later.
      typedef std::complex<real_Derived> Derived;
      typedef typename real_Derived::ancestor::cf_type real_cf_type;
      typedef std::complex<real_cf_type> cf_type;
      typedef typename real_Derived::ancestor::it_s_index real_it_s_index;
      typedef typename real_Derived::ancestor::term_type real_term_type;
      enum component
      {
        Real,
        Imag,
      };
// Build series from two components.
      void build_from_components(const real_Derived &p,
        const real_Derived &q)
      {
        p_assert(static_cast<Derived *>(this)->merge_args(p));
        if (!static_cast<Derived *>
          (this)->merge_args(q))
        {
          std::cout << "WARNING: constructing empty complex series because real and complex "
            "series passed to ctor are not argument compatible."<< std::endl;
          return;
        }
        typedef typename Derived::ancestor::term_type term_type;
// FIXME: hinted insertion here.
// Insert p (real_ part).
        {
          term_type term;
          real_it_s_index it=p.s_index().begin(), it_f=p.s_index().end();
          for (;
            it!=it_f;
            ++it)
          {
            term.s_cf()
              =cf_type(*it->g_cf(),real_cf_type(0));
            term.s_trig()=*it->g_trig();
            term.s_flavour()=it->g_flavour();
            static_cast<Derived *>(this)->insert(term);
          }
        }
// Insert q (imaginary part).
        {
          term_type term;
          real_it_s_index it=q.s_index().begin(), it_f=q.s_index().end();
          for (;
            it!=it_f;
            ++it)
          {
            term.s_cf()
              =cf_type(real_cf_type(0),*it->g_cf());
            term.s_trig()=*it->g_trig();
            term.s_flavour()=it->g_flavour();
            static_cast<Derived *>(this)->insert(term);
          }
        }
      }
      real_Derived get_comp(component cmp) const
      {
        typedef typename Derived::ancestor::it_s_index it_s_index;
        real_Derived retval;
        retval.merge_args(*static_cast<Derived const *>(this));
        retval.lin_args()=static_cast<Derived const *>(this)->lin_args();
        const it_s_index it_f=static_cast<Derived const *>(this)->s_index().end();
        real_it_s_index it_hint=retval.s_index().end();
        real_term_type term(real_cf_type(0),true);
        for (it_s_index it=static_cast<Derived const *>
          (this)->
          s_index().begin();
          it!=it_f;
          ++it)
        {
          term.s_cf()=get_cf_comp(*it->g_cf(),cmp);
          term.s_trig()=*it->g_trig();
          term.s_flavour()=it->g_flavour();
          it_hint=retval.insert(term,true,&it_hint);
        }
        return retval;
      }
// Extract component from complex pseries
      real_cf_type get_cf_comp(const cf_type &cf,component cmp) const
      {
        switch (cmp)
        {
          case Real:
            return cf.real();
          case Imag:
            return cf.imag();
        }
        std::cout << "WTF? enum is botched!" << std::endl;
        std::exit(1);
      }
    public:
/// Get real part.
      real_Derived real() const
      {
        return get_comp(Real);
      }
/// Get imaginary part.
      real_Derived imag() const
      {
        return get_comp(Imag);
      }
/// Absolute value.
      real_Derived abs() const
      {
        return (real()*real()+imag()*imag()).pow(.5);
      }
/// Complex conjugate.
      Derived conj() const
      {
        return Derived(real(),imag()*=-1);
      }
/// Make complex conjugate.
      void make_conj()
      {
        Derived tmp(real(),imag()*=-1);
        static_cast<Derived *>(this)->swap(tmp);
      }
  };
}
#endif
