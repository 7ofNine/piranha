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

#ifndef PIRANHA_BASE_PSERIES_DEF_H
#define PIRANHA_BASE_PSERIES_DEF_H

#include <boost/multi_index_container.hpp>
#include <boost/tuple/tuple.hpp>

#include "base_pseries_hooks.h"
#include "../bits/phase_list.h"
#include "../bits/psymbol.h"
#include "../bits/type_traits/eval_type.h"  // For eval_type.

namespace piranha
{
/// Template parameters for piranha::base_pseries.
#define __PIRANHA_BASE_PS_TP Cf,Trig,Term,I,Derived,Allocator
/// Template parameters for piranha::base_pseries (declaration form).
#define __PIRANHA_BASE_PS_TP_DECL class Cf, class Trig, template <class, class> class Term, template <class, class, \
  template <class, class> class> class I, class Derived, class Allocator

/// Base series class.
/**
 * Base class for the representation of Poisson series. It provides the methods for basic manipulation
 * of PS: I/O, management of terms, elementary maths, etc.
 *
 * This class should not be used directly, it should be inherited by a more specialized class.
 */
  template <class Cf, class Trig, template <class, class> class Term,
    template <class, class, template <class, class> class> class I, class Derived,
    class Allocator = std::allocator<char> >
    class base_pseries:base_pseries_hooks<base_pseries<__PIRANHA_BASE_PS_TP> >
  {
    public:
/// Alias for self.
      typedef base_pseries self;
/// Alias for coefficient type.
      typedef Cf cf_type;
/// Alias for trigonometric type.
      typedef Trig trig_type;
/// Alias for term type.
      typedef Term<cf_type,trig_type> term_type;
      typedef I<cf_type,trig_type,Term> index_type;
      typedef Derived derived_type;
      typedef Allocator allocator_type;
      typedef typename allocator_type::template rebind<term_type>::other term_allocator_type;
/// Alias for the evaluation type.
/**
 * @see base_pseries::t_eval.
 */
      typedef typename eval_type<cf_type>::type eval_type;
/// Alias for terms set type.
      typedef typename boost::multi_index_container <term_type,typename index_type::type,allocator_type>
        series_set_type;
/// Alias for sorted index.
      typedef typename series_set_type::template nth_index<0>::type sorted_index;
/// Alias for the iterator on sorted index.
      typedef typename sorted_index::const_iterator it_s_index;
/// Alias for the reversed iterator on sorted index.
      typedef typename sorted_index::const_reverse_iterator r_it_s_index;
/// Alias for hashed index.
      typedef typename series_set_type::template nth_index<1>::type hashed_index;
/// Alias for the iterator on hashed index.
      typedef typename hashed_index::const_iterator it_h_index;
/// Standard iterator.
/**
 * Standard iterator, so that piranha::base_pseries (and its children) can be used as STL
 * containers.
 * It can be defined here since it refers to base_pseries::series_set_type, which is the same object
 * in base_pseries and derived classes.
 * It also enables python iterators in pyranha.
 * @see base_pseries::begin().
 * @see base_pseries::end().
 */
      typedef it_s_index iterator;
      typedef it_s_index const_iterator;
// Ctors
      explicit base_pseries();
      explicit base_pseries(int);
      explicit base_pseries(const double &);
      base_pseries(const Derived &);
      explicit base_pseries(const std::string &);
      explicit base_pseries(const cf_type &, const Derived &);
      explicit base_pseries(const psymbol &, psymbol::type);
      ~base_pseries();
      Derived copy() const;
// Getters
      size_t length() const;
      size_t trig_width() const;
      size_t cf_width() const;
      const vector_int16 &lin_args() const;
      const vector_psym_p &cf_s_vec() const;
      const vector_psym_p &trig_s_vec() const;
      const series_set_type *g_series_set() const;
      const sorted_index &g_s_index() const;
      const hashed_index &g_h_index() const;
      int cf_arg_index(const std::string &) const;
      int trig_arg_index(const std::string &) const;
      size_t address() const;
      it_s_index begin() const;
      it_s_index end() const;
// Basic manipulation
      it_s_index insert(const term_type &term, bool sign = true, const it_s_index *it_hint = 0)
      {
        return insert_<true>(term,sign,it_hint);
      }
      it_s_index insert_no_sign_check(const term_type &term, bool sign = true, const it_s_index *it_hint = 0)
      {
        return insert_<false>(term,sign,it_hint);
      }
      template <class Cf2>
        it_s_index insert(const Term<Cf2, trig_type> &term,
        bool sign = true, const it_s_index *it_hint = 0)
      {
        return insert_<true>(term_type(term),sign,it_hint);
      }
      template <class Cf2>
        it_s_index insert_no_sign_check(const Term<Cf2, trig_type> &term,
        bool sign = true, const it_s_index *it_hint = 0)
      {
        return insert_<false>(term_type(term),sign,it_hint);
      }
      void term_erase(const it_h_index &);
      void term_erase(const it_s_index &);
      Derived &swap(Derived &);
      void cumulative_crop(const double &);
      void crop(const double &);
      void crop(const it_s_index &);
      void insert_phases(const phase_list &);
      void spectral_cutoff(const double &, const double &);
      template <class Derived2>
        bool merge_args(const Derived2 &);
/// Add coefficient argument.
/**
 * @param[in] psym, piranha::psymbol to be added.
 */
      void add_cf_arg(const psymbol &psym)
      {
        append_cf_args(vector_psym_p(1,psymbol_manager::get_pointer(psym)));
      }
/// Add trigonometric argument.
/**
 * @param[in] psym, piranha::psymbol to be added.
 */
      void add_trig_arg(const psymbol &psym)
      {
        append_trig_args(vector_psym_p(1,psymbol_manager::get_pointer(psym)));
      }
// I/O
      void print_plain(std::ostream &out_stream = std::cout, int limit = -1) const;
      void print_latex(std::ostream &out_stream = std::cout, int limit = -1) const;
/// Print series to std::ostream.
/**
 * Print first "limit" terms. If limit is negative, print all terms. The output format is read
 * from the piranha::stream_manager class.
 */
      void print(std::ostream &out_stream, int limit=-1) const
      {
        switch (stream_manager::format())
        {
          case stream_manager::plain:
            print_plain(out_stream,limit);
            break;
          case stream_manager::latex:
            print_latex(out_stream,limit);
        }
      }
/// Print to screen the first n terms, including series' header.
      void put(int n) const
      {
        print(std::cout,n);
      }
      void put() const
      {
        put(-1);
      }
/// Print to screen the first "limit" terms, without series' header.
      void put_terms(int limit) const
      {
        switch (stream_manager::format())
        {
          case stream_manager::plain:
            print_terms_plain(std::cout,limit);
            break;
          case stream_manager::latex:
            print_terms_latex(std::cout,limit);
        }
      }
      void put_terms() const
      {
        put_terms(-1);
      }
      void save_to(const std::string &) const;
      void put_phases_freqs(int limit) const;
      void put_phases_freqs() const
      {
        put_phases_freqs(-1);
      }
// Probing.
/// Check whether a series is empty or not.
      bool empty() const
      {
        return g_series_set()->empty();
      }
      double trig_density() const
      {
        const iterator it_f=end();
        double retval=0;
        if (empty())
        {
          return retval;
        }
        size_t count=0;
        for (iterator it=begin();it!=it_f;++it)
        {
          retval+=it->g_trig()->density(*this);
          ++count;
        }
        return (retval/count);
      }
      size_t footprint() const;
      double g_norm() const;
      it_s_index discontinuity() const;
      eval_type t_eval_brute(const double &) const;
      eval_type t_eval(const double &) const;
      size_t trig_index(const std::string &) const;
      eval_type mean(const double &, const double &,
        const size_t &n = 1000) const;
      bool checkup() const;
      bool is_cf() const;
    protected:
// Setters.
      vector_int16 &lin_args();
      vector_psym_p &cf_s_vec();
      vector_psym_p &trig_s_vec();
      series_set_type *s_series_set();
      sorted_index &s_s_index();
      hashed_index &s_h_index();
// Low level manipulation.
      void add_phase_to_term(const double &, iterator, term_type &, base_pseries &) const;
      void add_phase_to_term(const double &, const term_type &, term_type &, base_pseries &) const;
      void append_cf_args(const vector_psym_p &);
      void append_trig_args(const vector_psym_p &);
      void prepend_cf_args(const vector_psym_p &);
      void prepend_trig_args(const vector_psym_p &);
      it_h_index find_term(const term_type &t) const
      {
        return g_h_index().find(*t.g_trig());
      }
      it_s_index term_insert_new(const term_type &, bool, const it_s_index *it_hint);
      void term_update(const it_h_index &, cf_type &);
      it_s_index ll_insert(const term_type &, bool, const it_s_index *);
      template <class Cf2>
        it_s_index ll_insert(const Term<Cf2, trig_type> &,
        bool, const it_s_index *);
      template <bool CheckTrigSign>
        it_s_index insert_(const term_type &, bool sign = true, const it_s_index *it_hint = 0);
// Low level probing.
      it_s_index sdp_cutoff(const double &, const double &) const;
// Low level maths.
      template <class Derived2>
        Derived &assign_series(const Derived2 &);
      template <class Derived2, bool Sign>
        void alg_sum_lin_args(const Derived2 &);
      template <class Derived2, bool Sign>
        Derived &merge_with_series(const Derived2 &);
      template <class T>
        Derived &add_generic(const T &);
      template <class T>
        Derived &mult_by_generic(const T &);
      template <class T>
        Derived &divide_by_generic(const T &);
      template <class Derived2>
        bool series_multiplication_preliminaries(const Derived2 &, Derived &);
      template <class Derived2>
        bool series_multiplication_optimize_for_cf(const Derived2 &);
      template <class Cf2, class LightTermPair>
        static void term_by_term_multiplication_trig(const term_type &, const Term<Cf2,trig_type> &,
        LightTermPair &, cf_type &);
    public:
// Mathematics.
// Assignment.
      Derived &assign(const Derived &);
// Addition.
      template <class Derived2>
        Derived &add(const Derived2 &ps2)
      {
        return merge_with_series<Derived2,true>(ps2);
      }
      template <class Derived2>
        Derived &subtract(const Derived2 &ps2)
      {
        return merge_with_series<Derived2,false>(ps2);
      }
// Multiplication.
      Derived &mult_by(int);
      Derived &mult_by(const double &x)
      {
        return mult_by_generic(x);
      }
      template <class Derived2>
        Derived &mult_by_self(const Derived2 &);
      Derived &divide_by(int n)
      {
        return divide_by_generic(n);
      }
      Derived &divide_by(const double &x)
      {
        return divide_by_generic(x);
      }
    private:
/// Generic builder.
      template <class T>
        void generic_builder(const T &x)
      {
        term_type term = term_type();
        *(term.s_cf()) = cf_type(x);
        insert(term);
      }
      template <class Derived2>
        bool args_different(const Derived2 &) const;
      template <class Derived2>
        bool is_compatible(const Derived2 &) const;
// Low-level I/O.
      void print_terms_plain(std::ostream &, int ) const;
      void print_terms_latex(std::ostream &, int ) const;
      void read_data_from_file(std::ifstream &, const std::string &);
      void load_from(const std::string&);
      void identify_sections(std::ifstream &, const std::string &);
      void read_cf_arg(std::ifstream &);
      void read_trig_arg(std::ifstream &);
      void read_lin_args(std::ifstream &);
      void read_terms(std::ifstream &, const std::string &);
// Functors.
//Name comparison functor for psymbol pointers. Used in merging of arguments.
      struct psym_p_cmp
      {
        psym_p_cmp() {}
        bool operator()(psym_p p1, psym_p p2) const
        {
          return (p1->name()<p2->name());
        }
      };
/// Functor to update the coefficient.
      struct modifier_update_cf
      {
        modifier_update_cf(cf_type &new_cf):new_cf_(&new_cf) {}
        ~modifier_update_cf() {}
        void operator()(term_type &term)
        {
          term.s_cf()->swap(*new_cf_);
        }
// NOTICE: evaluate the impact of using const & here, esp. when using gmp
        cf_type *new_cf_;
      };
// Data members.
      vector_int16    lin_args_;
      vector_psym_p   cf_s_vec_;
      vector_psym_p   trig_s_vec_;
      series_set_type private_series_set_;
      static term_allocator_type term_allocator;
  };

  template <__PIRANHA_BASE_PS_TP_DECL>
    typename base_pseries<__PIRANHA_BASE_PS_TP>::term_allocator_type
    base_pseries<__PIRANHA_BASE_PS_TP>::term_allocator;
}

#endif
