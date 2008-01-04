/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
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

/** @file base_pseries_def.h
 \brief Definition of base Poisson series class.

 Longer description here.
 */

#include <boost/multi_index_container.hpp>
#include <boost/tuple/tuple.hpp>

#include "base_pseries_hooks.h"
#include "base_pseries_ta_macros.h"
#include "../phase_list.h"
#include "../../psymbol.h"
#include "../../type_traits/eval_type.h"  // For eval_type.
namespace piranha
{
  /// Base series class.
  /**
   * Base class for the representation of Poisson series. It provides the methods for basic manipulation
   * of PS: I/O, management of terms, elementary maths, etc.
   *
   * This class should not be used directly, it should be inherited by a more specialized class.
   */
  template <__PIRANHA_BASE_PS_TP_DECL= std::allocator<char> >
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
    /// Alias for index type.
    typedef I<cf_type,trig_type,Term> index_type;
    /// Alias for derived series.
    typedef Derived derived_type;
    /// Alias for allocator type.
    typedef Allocator allocator_type;
    /// Alias for allocator rebinding to term_type.
    typedef typename allocator_type::template rebind<term_type>::other term_allocator_type;
    /// Alias for the tuple of arguments vectors.
    typedef boost::tuple<vector_psym_p,vector_psym_p> arguments_tuple_type;
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
    /// Const counterpart of base_pseries::iterator.
    typedef it_s_index const_iterator;
    // Ctors.
    explicit base_pseries();
    explicit base_pseries(int);
    explicit base_pseries(const double &);
    base_pseries(const Derived &);
    explicit base_pseries(const std::string &);
    explicit base_pseries(const cf_type &, const Derived &);
    explicit base_pseries(const psymbol &, psymbol::type);
    ~base_pseries();
    Derived copy() const;
    // Getters.
    size_t length() const;
    const vector_psym_p &cf_args() const;
    const vector_psym_p &trig_args() const;
    size_t trig_width() const;
    size_t cf_width() const;
    const vector_int16 &lin_args() const;
    const arguments_tuple_type &arguments() const;
    const series_set_type *g_series_set() const;
    const sorted_index &g_s_index() const;
    const hashed_index &g_h_index() const;
    std::pair<bool,size_t> cf_arg_index(const std::string &) const;
    std::pair<bool,size_t> trig_arg_index(const std::string &) const;
    size_t address() const;
    it_s_index begin() const;
    it_s_index end() const;
    // Public manipulation.
    void swap(Derived &);
    void cumulative_crop(const double &);
    void crop(const double &);
    // Public I/O.
    void print(std::ostream &, int limit=-1) const;
    void put(int) const;
    void put() const;
    void put_terms(int) const;
    void put_terms() const;
    void save_to(const std::string &) const;
    void put_phases_freqs(int limit) const;
    void put_phases_freqs() const;
    // Public probing.
    bool is_empty() const;
    size_t footprint() const;
    double g_norm() const;
    eval_type t_eval_brute(const double &) const;
    eval_type t_eval(const double &) const;
    eval_type mean(const double &, const double &, const size_t &n = 1000) const;
    bool checkup() const;
    bool is_cf() const;
    // Public mathematics.
    // Assignment.
    Derived &assign(int);
    Derived &assign(const double &);
    Derived &assign(const Derived &);
    // Addition and subtraction.
    Derived &add(int);
    Derived &add(const double &);
    Derived &add(const Derived &);
    Derived &subtract(int);
    Derived &subtract(const double &);
    Derived &subtract(const Derived &);
    // Multiplication.
    Derived &mult_by(int);
    Derived &mult_by(const double &);
    Derived &mult_by(const Derived &);
    // Division.
    Derived &divide_by(int);
    Derived &divide_by(const double &);
    protected:
    // Protected ctors.
    template <class T>
    void generic_builder(const T &);
    // Setters.
    vector_int16 &lin_args();
    vector_psym_p &cf_args();
    vector_psym_p &trig_args();
    arguments_tuple_type &arguments();
    series_set_type *s_series_set();
    sorted_index &s_s_index();
    hashed_index &s_h_index();
    // Protected manipulation.
    void crop(const it_s_index &);
    void spectral_cutoff(const double &, const double &);
    template <class Derived2>
    void merge_args(const Derived2 &);
    template <class Term2, bool, bool>
    it_s_index insert(const Term2 &, const it_s_index *it_hint = 0);
    template <class Term2>
    it_s_index insert_with_checks(const Term2 &, const it_s_index *it_hint = 0);
    void term_erase(const it_h_index &);
    void term_erase(const it_s_index &);
    void append_cf_arg(const psym_p);
    void append_trig_arg(const psym_p);
    // Protected probing.
    it_s_index sdp_cutoff(const double &, const double &) const;
    it_s_index discontinuity() const;
    double trig_density() const;
    // Protected maths.
    template <class T>
    Derived &assign_generic(const T &);
    template <class Derived2>
    Derived &assign_series(const Derived2 &);
    template <class Derived2>
    Derived &add_series(const Derived2 &);
    template <class Derived2>
    Derived &subtract_series(const Derived2 &);
    template <class T>
    Derived &add_generic(const T &);
    template <class T>
    Derived &subtract_generic(const T &);
    template <class T>
    Derived &mult_by_generic(const T &);
    template <class Derived2>
    Derived &mult_by_series(const Derived2 &);
    template <class T>
    Derived &divide_by_generic(const T &);
    template <class Derived2>
    bool series_multiplication_preliminaries(const Derived2 &, Derived &);
    template <class Derived2>
    bool series_multiplication_optimize_for_cf(const Derived2 &);
    template <class Cf2, class LightTermPair>
    static void term_by_term_multiplication_trig(const term_type &, const Term<Cf2,trig_type> &,
      LightTermPair &, cf_type &);
    private:
    // Private getters.
    template <int>
    std::pair<bool,size_t> arg_index(const std::string &) const;
    // Private probing.
    template <class Derived2>
    bool is_args_compatible(const Derived2 &) const;
    template <int>
    static std::string psymbol_descr();
    // Private manipulation.
    template <int>
    void append_arg(const psym_p);
    template <class Derived2>
    void merge_incompatible_args(const Derived2 &);
    it_h_index find_term(const term_type &) const;
    template <bool>
    it_s_index term_insert_new(const term_type &, const it_s_index *it_hint);
    void term_update(const it_h_index &, cf_type &);
    template <bool>
    it_s_index ll_insert(const term_type &, const it_s_index *);
    void add_phase_to_term(const double &, iterator, term_type &, base_pseries &) const;
    void add_phase_to_term(const double &, const term_type &, term_type &, base_pseries &) const;
    void insert_phases(const phase_list &);
    // Private maths.
    template <class Derived2, bool>
    void alg_sum_lin_args(const Derived2 &);
    template <class Derived2, bool>
    Derived &merge_with_series(const Derived2 &);
    // Private I/O.
    void print_plain(std::ostream &out_stream = std::cout, int limit = -1) const;
    void print_latex(std::ostream &out_stream = std::cout, int limit = -1) const;
    void print_terms_plain(std::ostream &, int) const;
    void print_terms_latex(std::ostream &, int) const;
    void read_data_from_file(std::ifstream &, const std::string &);
    void load_from(const std::string &);
    void identify_sections(std::ifstream &, const std::string &);
    template <psymbol::type>
    void read_arg(std::ifstream &);
    void read_cf_arg(std::ifstream &);
    void read_trig_arg(std::ifstream &);
    void read_lin_args(std::ifstream &);
    void read_terms(std::ifstream &, const std::string &);
    // Functors.
    // Functor to update the coefficient.
    struct modifier_update_cf
    {
      modifier_update_cf(cf_type &new_cf):new_cf_(&new_cf)
      {}
      ~modifier_update_cf()
      {}
      void operator()(term_type &term)
      {
        term.cf().swap(*new_cf_);
      }
      cf_type *new_cf_;
    };
    // Data members.
    vector_int16 lin_args_;
    arguments_tuple_type m_arguments;
    series_set_type private_series_set_;
    static term_allocator_type term_allocator;
  };

  template <__PIRANHA_BASE_PS_TP_DECL>
  typename base_pseries<__PIRANHA_BASE_PS_TP>::term_allocator_type
  base_pseries<__PIRANHA_BASE_PS_TP>::term_allocator;
}

#endif
