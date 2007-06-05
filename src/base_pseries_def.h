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
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/member.hpp>
#include <complex>

#include "phase_list.h"
#include "ps_term.h"
#include "ps_operators.h"

namespace piranha
  {
  /// Default indices for base_pseries.
  /**
   * This class specifies the default indices used in Poisson series: an hashed index for the identification
   * of terms and a sorted index to discard terms in multiplications. The class is to be used as the I
   * parameter in piranha::base_pseries classes.
   */
  template <class Cf, class Trig>
  struct default_ps_index
    {
      typedef boost::multi_index::indexed_by <
      boost::multi_index::ordered_unique <
      boost::multi_index::identity<ps_term<Cf, Trig> > > ,
      boost::multi_index::hashed_unique <
      boost::multi_index::composite_key <
      ps_term<Cf, Trig>,
      boost::multi_index::const_mem_fun < ps_term<Cf, Trig>, bool,
      &ps_term<Cf, Trig>::flavour > ,
      boost::multi_index::const_mem_fun < ps_term<Cf, Trig>, const Trig &,
      &ps_term<Cf, Trig>::trig_args >
      >
      >
      > type;
    };


  /// Base series class.
  /**
   * Base class for the representation of Poisson series (PS). It provides the methods for basic manipulation
   * of PS: I/O, management of terms, elementary maths, etc.
   *
   * This class should not be used directly, it should be inherited by a more specialized class.
   * A default specialized class, piranha::ps, exists.
   * @see piranha:ps, default specialized Poisson series class.
   */
  template <class Cf, class Trig, template <class, class> class I>
  class base_pseries
    {
    public:
      /// Alias for itself.
      typedef base_pseries self;
      /// Alias for coefficient type.
      typedef Cf cf_type;
      /// Alias for the evaluation type.
      /**
       * @see base_pseries::t_eval.
       */
      typedef typename cf_type::eval_type eval_type;
      /// Alias for trigonometric type.
      typedef Trig trig_type;
      typedef I<cf_type, trig_type> index_type;
      /// Alias for term type.
      typedef ps_term<cf_type, trig_type> term_type;
      typedef typename boost::multi_index_container < term_type,
      typename index_type::type > set_type;
      /// Alias for the sorted index.
      typedef typename set_type::template nth_index<0>
      ::type sorted_index;
      /// Alias for the iterator on sorted index.
      typedef typename sorted_index::iterator it_s_index;
      /// Alias for the reversed iterator on sorted index.
      typedef typename sorted_index::const_reverse_iterator r_it_s_index;
      /// Alias for the hashed index.
      typedef typename set_type::template nth_index<1>
      ::type hashed_index;
      /// Alias for the iterator on hashed index.
      typedef typename hashed_index::iterator it_h_index;
      /// Standard iterator.
      /**
       * Standard iterator, so that piranha::base_pseries (and its children) can be used as STL
       * containers.
       * It can be defined here since it refers to base_pseries::set_type, which is the same object
       * in base_pseries and derived classes.
       * It also enables python iterators in pyranha.
       * @see base_pseries::begin().
       * @see base_pseries::end().
       * @see base_pseries::const_iterator.
       */
      typedef it_s_index iterator;
      /// Standard constant iterator.
      /**
       * Constant counterpart of base_pseries::iterator.
       * @see base_pseries::iterator.
       */
      typedef typename sorted_index::const_iterator const_iterator;
      // Ctors
      base_pseries();
      base_pseries(const base_pseries &);
      explicit base_pseries(const std::string &);
      explicit base_pseries(const cf_type &);
      explicit base_pseries(const psymbol &);
      /// Destructor.
      ~base_pseries()
      {}
      // Getters
      /// Return series length.
      size_t length() const
        {
          return set_.size();
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
      /// Return a reference to the set of terms.
      const set_type &set
      () const
        {
          return set_;
        }
      psymbol trig_arg(const size_t &) const;
      psymbol trig_arg(const std::string &) const;
      /// Return a const reference to the sorted index.
      /**
       * @see base_pseries::sorted_index.
       */
      const sorted_index &s_index() const
        {
          return set_.template get
                   <0>();
        }
      /// Return a const reference to the hashed index.
      /**
       * @see base_pseries::hashed_index.
       */
      const hashed_index &h_index() const
        {
          return set_.template get
                   <1>();
        }
      // Setters for indices.
      sorted_index &s_index()
      {
        return set_.template get
                 <0>();
      }
      hashed_index &h_index()
      {
        return set_.template get
                 <1>();
      }
      /// Return a numerical value corresponding to the memory address of the series.
      size_t address()
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
          return s_index().begin();
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
          return s_index().end();
        }
      // Basic manipulation
      it_s_index insert(const term_type &, bool sign = true, const it_s_index *it_hint = 0);
      template <class Cf2, class Trig2>
      it_s_index insert(const ps_term<Cf2, Trig2> &,
                        bool sign = true, const it_s_index *it_hint = 0);
      void swap(base_pseries &);
      void cumulative_crop(const double &);
      void crop(const double &);
      void crop(const it_s_index &);
      void insert_phases(const phase_list &);
      void spectral_cutoff(const double &, const double &);
      template <class Cf2, class Trig2, template <class, class> class I2>
      bool merge_args(const
                      base_pseries<Cf2, Trig2, I2> &);
      void set_flavour(bool);
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
      // Maths.
      __PS_OPERATORS(base_pseries, base_pseries);
      base_pseries &operator*=(const cf_type &x)
      {
        generic_mult(x);
        return *this;
      }
      // Probing.
      /// Check whether a series is empty or not.
      bool empty() const
        {
          return set_.empty();
        }
      /// Return series' norm.
      /**
       * The norm is kept up-to-date during term insertions. It is calculated using a norm() method
       * provided by the coefficient class.
       */
      double norm() const
        {
          return norm_;
        }
      size_t footprint() const;
      double calc_norm() const;
      it_s_index discontinuity() const;
      eval_type t_eval(const double &) const;
      size_t trig_index(const std::string &) const;
      eval_type mean(const double &, const double &,
                     const size_t &n = 1000) const;
      bool checkup() const;
      // FIXME: temporarily here.
    protected:
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
      void upgrade_norm(const double &);
      void downgrade_norm(const double &);
      it_s_index term_insert_new(const term_type &, bool, const it_s_index *it_hint);
      void term_erase(const it_h_index &);
      void term_erase(const it_s_index &);
      void term_update(const it_h_index &, const cf_type &);
      void term_update(const it_s_index &, const cf_type &);
      it_s_index ll_insert(const term_type &, bool, const it_s_index *);
      template <class Cf2, class Trig2>
      it_s_index ll_insert(const ps_term<Cf2, Trig2> &,
                           bool, const it_s_index *);
      // Low level probing.
      it_s_index sdp_cutoff(const double &, const double &) const;
      template <class Cf2, class Trig2, template <class, class> class I2>
      bool args_compatible(const
                              base_pseries<Cf2, Trig2, I2> &) const;
      struct psym_p_cmp
      {
        psym_p_cmp()
        {}
        bool operator()(psym_p p1, psym_p p2) const
          {
            return (p1->name()<p2->name());
          }
      };
      template <class Cf2, class Trig2, template <class, class> class I2>
      bool args_different(const
                              base_pseries<Cf2, Trig2, I2> &) const;
      // Low level maths.
      void basic_assignment(const base_pseries &);
      template <class Cf2, class Trig2, template <class, class> class I2>
      void alg_sum_lin_args(const base_pseries<Cf2, Trig2, I2> &, bool);
      template <class Cf2, class Trig2, template <class, class> class I2>
      void merge_with(const base_pseries<Cf2, Trig2, I2>&, bool sign = true);
      template <class T>
      void generic_merge(const T &);
      template <class Cf2, class Trig2, template <class, class> class I2>
      void mult_terms(const base_pseries<Cf2, Trig2, I2> &, base_pseries &, const double &);
      template <class T>
      void basic_ps_mult(const T &);
      template <class T>
      void generic_mult(const T &);
      void mult_by_int(int);
      unsigned int pow_limit(const double &, const double &) const;
      void basic_pow(const double &);
      // Data members.
    protected:
      double          norm_;
      vector_mult_t   lin_args_;
      vector_psym_p   cf_s_vec_;
      vector_psym_p   trig_s_vec_;
      set_type        set_;
    };


  // Default ctor
#define __base_pseries_init_list norm_(0),lin_args_(),cf_s_vec_(),trig_s_vec_(),set_()

  /// Default constructor.
  /**
   * Constructs an empty series.
   */
  template <class Cf, class Trig, template <class, class> class I>
  inline base_pseries<Cf, Trig, I>::base_pseries(): __base_pseries_init_list {}


  /// Copy constructor.
  /**
   * Constructs a series from another one.
   */
  template <class Cf, class Trig, template <class, class> class I>
  inline base_pseries<Cf, Trig, I>::base_pseries(const base_pseries &ps):
      norm_(ps.norm_), lin_args_(ps.lin_args_), cf_s_vec_(ps.cf_s_vec_),
      trig_s_vec_(ps.trig_s_vec_), set_(ps.set_)
  {
    std::cout << "Copy ctor" << std::endl;
  }


  /// Constructor from filename.
  /**
   * Read a series from file.
   */
  template <class Cf, class Trig, template <class, class> class I>
  inline base_pseries<Cf, Trig, I>::base_pseries(const std::string &fn): __base_pseries_init_list
  {
    load_from(fn);
  }



  /// Constructor from coefficient.
  /**
   * Constructs a series consisting of a single cosine term with zero trigonometric arguments
   * and provided coefficient.
   * @see base_pseries::cf_type.
   */
  template <class Cf, class Trig, template <class, class> class I>
  inline base_pseries<Cf, Trig, I>::base_pseries(const cf_type &c): __base_pseries_init_list
  {
    term_type term(c);
    insert(term);
  }


  /// Constructor from piranha::psymbol.
  template <class Cf, class Trig, template <class, class> class I>
  inline base_pseries<Cf, Trig, I>::base_pseries(const psymbol &psym): __base_pseries_init_list
  {
    append_cf_args(vector_psym_p(1,psymbol_manager::get_pointer(psym)));
    cf_type c(psym);
    term_type term(c);
    insert(term);
  }

#undef __base_pseries_init_list


  /// Return a copy of the trigonometric argument corresponding to a certain position n.
  template <class Cf, class Trig, template <class, class> class I>
  inline psymbol base_pseries<Cf, Trig, I>::trig_arg(const size_t &n) const
    {
      if (n >= trig_width())
        {
          std::cout << "WARNING: invalid arg index requested, returning null argument." << std::endl;
          return psymbol();
        }
      return *trig_s_vec_[n];
    }


  /// Return a copy of the trigonometric argument by name s.
  template <class Cf, class Trig, template <class, class> class I>
  inline psymbol base_pseries<Cf, Trig, I>::trig_arg(const std::string &s) const
    {
      return trig_arg(trig_index(s));
    }
}
#endif
