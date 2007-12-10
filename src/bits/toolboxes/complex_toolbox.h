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
/// Complex toolbox.
  template <class real_Derived>
    class complex_toolbox
  {
/// Alias for complex type.
      typedef std::complex<real_Derived> Derived;
    public:
      explicit complex_toolbox() {}
// Complex specific constructors.
      explicit complex_toolbox(const std::complex<int> &n)
      {
        static_cast<Derived *>(this)->generic_builder(n);
      }
      explicit complex_toolbox(const std::complex<double> &x)
      {
        static_cast<Derived *>(this)->generic_builder(x);
      }
      explicit complex_toolbox(int r, int i)
      {
        static_cast<Derived *>(this)->generic_builder(std::complex<int>(r,i));
      }
      explicit complex_toolbox(const double &r, const double &i)
      {
        static_cast<Derived *>(this)->generic_builder(std::complex<double>(r,i));
      }
      ~complex_toolbox() {}
/// Get real part.
      real_Derived real() const
      {
        return get_comp<Real>();
      }
/// Get imaginary part.
      real_Derived imag() const
      {
        return get_comp<Imag>();
      }
/// Absolute value.
// TODO:place this into pow toolbox, complex counterpart?
      real_Derived abs() const
      {
        return (real()*real()+imag()*imag()).pow(.5);
      }
/// Complex conjugate.
      Derived conj() const
      {
        return Derived(real(),imag()*=-1);
      }
// Maths interoperability with complex pods and real counterpart.
/// Assign complex int.
      Derived &complex_assign(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->assign_generic(n);
      }
/// Assign complex double.
      Derived &complex_assign(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->assign_generic(x);
      }
/// Assign real counterpart.
      Derived &complex_assign(const real_Derived &r)
      {
        return static_cast<Derived *>(this)->assign_series(r);
      }
/// Add complex int.
      Derived &complex_add(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->add_generic(n);
      }
/// Add complex double.
      Derived &complex_add(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->add_generic(x);
      }
/// Add real counterpart.
      Derived &complex_add(const real_Derived &r)
      {
        return static_cast<Derived *>(this)->add_series(r);
      }
/// Subtract complex int.
      Derived &complex_subtract(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->subtract_generic(n);
      }
/// Subtract complex double.
      Derived &complex_subtract(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->subtract_generic(x);
      }
/// Subtract real counterpart.
      Derived &complex_subtract(const real_Derived &r)
      {
        return static_cast<Derived *>(this)->subtract_series(r);
      }
/// Multiply by complex int.
      Derived &complex_mult_by(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->mult_by_generic(n);
      }
/// Multiply by complex double.
      Derived &complex_mult_by(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->mult_by_generic(x);
      }
/// Multiply by real counterpart.
      Derived &complex_mult_by(const real_Derived &r)
      {
        return static_cast<Derived *>(this)->mult_by_series(r);
      }
/// Divide by complex int.
      Derived &complex_divide_by(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->divide_by_generic(n);
      }
/// Divide by complex double.
      Derived &complex_divide_by(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->divide_by_generic(x);
      }
    protected:
// NOTICE: typedefs regarding Derived type cannot be placed here because when the compiler
// parses this part it does not know enough about Derived yet. They are ok in the body of
// the methods tough, because instantiation happens later.
      typedef typename real_Derived::ancestor::cf_type real_cf_type;
      typedef std::complex<real_cf_type> cf_type;
      typedef typename real_Derived::ancestor::it_s_index real_it_s_index;
      typedef typename real_Derived::ancestor::term_type real_term_type;
// TODO: what's with the extra comma at the end?
      enum component
      {
        Real,
        Imag,
      };
      template <component Cmp>
        void insert_component(const real_Derived &comp)
      {
        typedef typename Derived::ancestor::term_type term_type;
        term_type term;
        const real_it_s_index it_f=comp.g_s_index().end();
// TODO: hinted insertion here.
        for (real_it_s_index it=comp.g_s_index().begin();it!=it_f;++it)
        {
          *term.s_cf()=build_cf_from_comp<Cmp>(*it->g_cf());
          *term.s_trig()=*it->g_trig();
          term.s_flavour()=it->g_flavour();
          static_cast<Derived *>(this)->insert(term);
        }
      }
// Build series from real and imaginary components.
      void build_from_components(const real_Derived &p, const real_Derived &q)
      {
        action_assert(static_cast<Derived *>(this)->merge_args(p));
        if (!static_cast<Derived *>(this)->merge_args(q))
        {
          std::cout << "WARNING: constructing empty complex series because real and complex "
            "series passed to ctor are not argument compatible."<< std::endl;
          return;
        }
        insert_component<Real>(p);
        insert_component<Imag>(q);
      }
      template <component Cmp>
        real_Derived get_comp() const
      {
        typedef typename Derived::ancestor::it_s_index it_s_index;
        real_Derived retval;
        retval.merge_args(*static_cast<Derived const *>(this));
        retval.lin_args()=static_cast<Derived const *>(this)->lin_args();
        const it_s_index it_f=static_cast<Derived const *>(this)->g_s_index().end();
        real_it_s_index it_hint=retval.g_s_index().end();
        real_term_type term(real_cf_type(0));
        for (it_s_index it=static_cast<Derived const *>(this)->g_s_index().begin();it!=it_f;++it)
        {
          *term.s_cf()=get_cf_comp<Cmp>(*it->g_cf());
          *term.s_trig()=*it->g_trig();
          term.s_flavour()=it->g_flavour();
          it_hint=retval.insert(term,&it_hint);
        }
        return retval;
      }
    private:
// Extract component from complex pseries
      template <component Cmp>
        real_cf_type get_cf_comp(const cf_type &cf) const
      {
        switch (Cmp)
        {
          case Real:
            return cf.real();
          default:
            return cf.imag();
        }
      }
      template <component Cmp>
        cf_type build_cf_from_comp(const real_cf_type &cf) const
      {
        switch (Cmp)
        {
          case Real:
            return cf_type(cf,real_cf_type(0));
          default:
            return cf_type(real_cf_type(0),cf);
        }
      }
  };
}
#endif
