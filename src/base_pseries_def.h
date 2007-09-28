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
#include <ext/pool_allocator.h>

#include "base_pseries_hooks.h"
#include "phase_list.h"
#include "psymbol.h"

namespace piranha
{
/// Base series class.
/**
 * Base class for the representation of Poisson series (PS). It provides the methods for basic manipulation
 * of PS: I/O, management of terms, elementary maths, etc.
 *
 * This class should not be used directly, it should be inherited by a more specialized class.
 * A default specialized class, piranha::ps, exists.
 * @see piranha:ps, default specialized Poisson series class.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    class base_pseries:base_pseries_hooks<base_pseries<Cf,Trig,Term,I,Derived> >
  {
    public:
/// Alias for self.
      typedef base_pseries self;
/// Alias for coefficient type.
      typedef Cf cf_type;
/// Alias for trigonometric type.
      typedef Trig trig_type;
/// Alias for term type.
      typedef Term<cf_type, trig_type> term_type;
      typedef I<cf_type, trig_type, Term> index_type;
      typedef Derived derived_type;
/// Alias for the evaluation type.
/**
 * @see base_pseries::t_eval.
 */
      typedef typename cf_type::eval_type eval_type;
      typedef typename boost::multi_index_container < term_type,
        typename index_type::type, __gnu_cxx::__pool_alloc<term_type> > series_set_type;
/// Alias for the sorted index.
      typedef typename series_set_type::template nth_index<0>
        ::type sorted_index;
/// Alias for the iterator on sorted index.
      typedef typename sorted_index::const_iterator it_s_index;
/// Alias for the reversed iterator on sorted index.
      typedef typename sorted_index::const_reverse_iterator r_it_s_index;
/// Alias for the hashed index.
      typedef typename series_set_type::template nth_index<1>
        ::type hashed_index;
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
      base_pseries();
      base_pseries(const Derived &);
      explicit base_pseries(const std::string &);
      explicit base_pseries(const cf_type &, const Derived &);
      explicit base_pseries(const psymbol &, psymbol::type);
/// Destructor.
      ~base_pseries()
        {}
/// Copy.
      Derived copy() const
      {
        return Derived(*static_cast<Derived const *>(this));
      }
// Getters
/// Return series length.
      size_t length() const
      {
        return g_series_set()->size();
      }
/// Return trigonometric width.
      size_t trig_width() const
      {
//p_assert(trig_s_vec_.size()==lin_args_.size());
        return trig_s_vec_.size();
      }
/// Return coefficient width.
      size_t cf_width() const
      {
        return cf_s_vec_.size();
      }
/// Return reference to the vector of linear arguments.
      vector_mult_t &lin_args()
      {
        return lin_args_;
      }
      const vector_mult_t &lin_args() const
      {
        return lin_args_;
      }
/// Return reference to the vector of coefficient symbols.
      vector_psym_p &cf_s_vec()
      {
        return cf_s_vec_;
      }
      const vector_psym_p &cf_s_vec() const
      {
        return cf_s_vec_;
      }
/// Return a reference to the vector of trigonometric symbols.
      vector_psym_p &trig_s_vec()
      {
        return trig_s_vec_;
      }
      const vector_psym_p &trig_s_vec() const
      {
        return trig_s_vec_;
      }
/// Return a const reference to the set of terms.
      const series_set_type *g_series_set() const
      {
        return &private_series_set_;
      }
/// Return a reference to the set of terms.
      series_set_type *s_series_set()
      {
        return &private_series_set_;
      }
/// Return a const reference to the sorted index.
/**
 * @see base_pseries::sorted_index.
 */
      const sorted_index &g_s_index() const
      {
        return g_series_set()->template get
          <0>();
      }
/// Return a const reference to the hashed index.
/**
 * @see base_pseries::hashed_index.
 */
      const hashed_index &g_h_index() const
      {
        return g_series_set()->template get
          <1>();
      }
/// Return index of coefficient argument "name".
      int cf_arg_index(const std::string &name) const
      {
        int retval=-1;
        const size_t w=cf_width();
        for (size_t i=0;i<w;++i)
        {
          if (cf_s_vec_[i]->name()==name)
          {
            retval=(int)i;
            break;
          }
        }
        return retval;
      }
/// Return index of trigonometric argument "name".
      int trig_arg_index(const std::string &name) const
      {
        int retval=-1;
        const size_t w=trig_width();
        for (size_t i=0;i<w;++i)
        {
          if (trig_s_vec_[i]->name()==name)
          {
            retval=(int)i;
            break;
          }
        }
        return retval;
      }
// Setters for indices.
      sorted_index &s_s_index()
      {
        return s_series_set()->template get
          <0>();
      }
      hashed_index &s_h_index()
      {
        return s_series_set()->template get
          <1>();
      }
/// Return a numerical value corresponding to the memory address of the series.
      size_t address() const
      {
        return (size_t)(void *)this;
      }
/// Begin of series.
/**
 * Returns an iterator pointing to the first term of the series. This allows to mimic the behaviour
 * of an STL container.
 * @see base_pseries::end().
 * @see base_pseries::iterator.
 */
      it_s_index begin() const
      {
        return g_s_index().begin();
      }
/// End of series.
/**
 * Returns an iterator pointing to the last term of the series. This allows to mimic the behaviour
 * of an STL container.
 * @see base_pseries::begin().
 * @see base_pseries::iterator.
 */
      it_s_index end() const
      {
        return g_s_index().end();
      }
// Basic manipulation
      it_s_index insert(const term_type &, bool sign = true, const it_s_index *it_hint = 0);
      template <class Cf2>
        it_s_index insert(const Term<Cf2, trig_type> &,
        bool sign = true, const it_s_index *it_hint = 0);
      void term_erase(const it_h_index &);
      void term_erase(const it_s_index &);
      void swap(Derived &);
      void cumulative_crop(const double &);
      void crop(const double &);
      void crop(const it_s_index &);
      void insert_phases(const phase_list &);
      void spectral_cutoff(const double &, const double &);
      template <class Derived2>
        bool merge_args(const Derived2 &);
      void set_flavour(bool);
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
/// Print to screen the first "limit" terms, including series' header.
      void put(int limit) const
      {
        print(std::cout, limit);
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
// NOTICE: temporarily here.
      template <class Derived2>
        void series_multiplication(const Derived2 &);
      template <class Derived2>
        void generic_series_assignment(const Derived2 &);
    protected:
/// Generic builder.
      template <class T>
        void generic_builder(const T &x)
      {
        term_type term(x);
        insert(term);
      }
// Low level I/O.
      void read_data_from_file(std::ifstream &, const std::string &);
      void load_from(const std::string&);
      void identify_sections(std::ifstream &, const std::string &);
      void read_cf_arg(std::ifstream &);
      void read_trig_arg(std::ifstream &);
      void read_lin_args(std::ifstream &);
      void read_terms(std::ifstream &, const std::string &);
      void print_terms_plain(std::ostream &, int ) const;
      void print_terms_latex(std::ostream &, int ) const;
// Low level manipulation.
      void add_phase_to_term(const double &, iterator, term_type &, base_pseries &) const;
      void add_phase_to_term(const double &, const term_type &, term_type &, base_pseries &) const;
      void append_cf_args(const vector_psym_p &);
      void append_trig_args(const vector_psym_p &);
      void prepend_cf_args(const vector_psym_p &);
      void prepend_trig_args(const vector_psym_p &);
      it_s_index term_insert_new(const term_type &, bool, const it_s_index *it_hint);
      void term_update(const it_h_index &, cf_type &);
      it_s_index ll_insert(const term_type &, bool, const it_s_index *);
      template <class Cf2>
        it_s_index ll_insert(const Term<Cf2, trig_type> &,
        bool, const it_s_index *);
// Low level probing.
      it_s_index sdp_cutoff(const double &, const double &) const;
// Low level maths.
      void basic_assignment(const base_pseries &);
      template <class Derived2>
        void alg_sum_lin_args(const Derived2 &, bool);
    public:
      it_h_index find_term(const term_type &t) const
      {
        return g_h_index().find(*t.g_trig());
      }
      template <class Derived2>
        void merge_with(const Derived2 &, bool sign = true);
      template <class T>
        void generic_merge(const T &);
      template <class T>
        void generic_multiplication(const T &);
      template <class T>
        void generic_division(const T &);
      void mult_by_int(int);
      template <class Derived2>
        bool series_multiplication_preliminaries(const Derived2 &, Derived &);
      template <class Derived2>
        bool series_multiplication_optimize_for_cf(const Derived2 &);
      template <class Cf2, class V>
        static void term_by_term_multiplication_trig(const term_type &, const Term<Cf2,trig_type> &,
        V &, cf_type &);
    private:
/// Name comparison functor for psymbol pointers.
/**
 * Used in merging of arguments.
 */
      struct psym_p_cmp
      {
        psym_p_cmp()
          {}
        bool operator()(psym_p p1, psym_p p2) const
        {
          return (p1->name()<p2->name());
        }
      };
      template <class Cf2, class Derived2>
        bool args_different(const
        base_pseries<Cf2, trig_type, Term, I, Derived2> &) const;
      template <class Derived2>
        bool args_compatible(const Derived2 &) const;
    private:
/// Functor to update the coefficient.
// TODO: use swap here?
      struct modifier_update_cf
      {
        modifier_update_cf(cf_type &new_cf):new_cf_(&new_cf)
          {}
        ~modifier_update_cf()
          {}
        void operator()(term_type &term)
        {
          term.s_cf()->swap(*new_cf_);
        }
// NOTICE: evaluate the impact of using const & here, esp. when using gmp
        cf_type *new_cf_;
      };
// Data members.
    protected:
      vector_mult_t   lin_args_;
      vector_psym_p   cf_s_vec_;
      vector_psym_p   trig_s_vec_;
      series_set_type private_series_set_;
  };

// Default ctor
#define __base_pseries_init_list lin_args_(),cf_s_vec_(),trig_s_vec_(),private_series_set_()

/// Default constructor.
/**
 * Constructs an empty series.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline base_pseries<Cf, Trig, Term, I, Derived>::base_pseries():__base_pseries_init_list {}

/// Copy constructor.
/**
 * Constructs a series from another one.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline base_pseries<Cf, Trig, Term, I, Derived>::base_pseries(const Derived &ps):
  lin_args_(ps.lin_args_),cf_s_vec_(ps.cf_s_vec_),trig_s_vec_(ps.trig_s_vec_),private_series_set_(*ps.g_series_set())
  {
    std::cout << "Copy ctor" << std::endl;
  }

/// Constructor from filename.
/**
 * Read a series from file.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline base_pseries<Cf, Trig, Term, I, Derived>::base_pseries(const std::string &fn):__base_pseries_init_list
  {
    load_from(fn);
  }

/// Constructor from coefficient.
/**
 * Constructs a series consisting of a single cosine term with zero trigonometric arguments
 * and provided coefficient.
 * @see base_pseries::cf_type.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline base_pseries<Cf, Trig, Term, I, Derived>::base_pseries(const cf_type &c, const Derived &model)
    :__base_pseries_init_list
  {
    if (!merge_args(model))
    {
      std::cout << "Warning: incompatbile arguments in ctor from cf_type." << std::endl;
      return;
    }
    if (c.larger(cf_width()))
    {
      std::cout << "Warning: too many arguments in ctor from coefficient." << std::endl;
      return;
    }
    generic_builder(c);
  }

/// Constructor from piranha::psymbol.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline base_pseries<Cf, Trig, Term, I, Derived>::base_pseries(const psymbol &psym, psymbol::type ptype)
    :__base_pseries_init_list
  {
    if (ptype==psymbol::cf)
    {
// When building to cf create a coefficient from the symvol.
      append_cf_args(vector_psym_p(1,psymbol_manager::get_pointer(psym)));
      cf_type c(psym);
      term_type term(c);
      insert(term);
    }
    else if (ptype==psymbol::trig)
    {
// When building to trig assign argument in lin_args.
      append_trig_args(vector_psym_p(1,psymbol_manager::get_pointer(psym)));
      lin_args_[0]=1;
    }
    else
    {
      p_assert(false);
    }
  }

#undef __base_pseries_init_list

}
#endif
