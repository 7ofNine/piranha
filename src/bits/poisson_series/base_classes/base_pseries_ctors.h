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

#ifndef PIRANHA_BASE_PSERIES_CTORS_H
#define PIRANHA_BASE_PSERIES_CTORS_H

#include "../../p_exceptions.h"
#include "base_pseries_ta_macros.h"

namespace piranha
{
/// Initializer list for base_pseries constructors.
#define __base_pseries_init_list lin_args_(),m_arguments(),private_series_set_()

/// Default constructor.
/**
 * Constructs a null series: empty with zero arguments.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries():__base_pseries_init_list {}

/// Constructor from int.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries(int n):__base_pseries_init_list
  {
    generic_builder(n);
  }

/// Constructor from double.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries(const double &x):__base_pseries_init_list
  {
    generic_builder(x);
  }

/// Copy constructor.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries(const Derived &ps):lin_args_(ps.lin_args_),
    m_arguments(ps.m_arguments),private_series_set_(*ps.g_series_set())
  {
    std::cout << "Copy ctor" << std::endl;
  }

/// Constructor from file.
/**
 * Read a series from file.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries(const std::string &fn):__base_pseries_init_list
  {
    load_from(fn);
  }

/// Constructor from coefficient.
/**
 * Constructs a series consisting of a single cosine term with zero trigonometric arguments,
 * provided coefficient and same set of arguments as model.
 * @see base_pseries::cf_type.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries(const cf_type &c, const Derived &model):__base_pseries_init_list
  {
    try
    {
      merge_args(model);
      generic_builder(c);
    }
    catch (exceptions::add_arguments &e)
    {
      std::cout << "Too many coefficient arguments in ctor from coefficient, building empty series." << std::endl;
    }
  }

/// Constructor from piranha::psymbol.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::base_pseries(const psymbol &psym, psymbol::type ptype):__base_pseries_init_list
  {
    switch (ptype)
    {
      case (psymbol::cf):
        {
// When building as cf create a coefficient from the symbol.
          append_cf_args(vector_psym_p(1,psymbol_manager::get_pointer(psym).second));
          cf_type c(psym);
          term_type term(c,trig_type());
          insert_check_positive(term);
        }
        break;
      case (psymbol::trig):
        {
// When building as trig assign argument in lin_args.
          append_trig_args(vector_psym_p(1,psymbol_manager::get_pointer(psym).second));
          lin_args_[0]=1;
        }
        break;
      default:
        p_assert(false);
        ;
    }
  }

/// Destructor.
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline base_pseries<__PIRANHA_BASE_PS_TP>::~base_pseries() {}

#undef __base_pseries_init_list

/// Copy.
/**
 * Returns a copy of this.
 */
  template <__PIRANHA_BASE_PS_TP_DECL>
    inline Derived base_pseries<__PIRANHA_BASE_PS_TP>::copy() const
  {
    return Derived(*static_cast<Derived const *>(this));
  }

/// Generic builder.
  template <__PIRANHA_BASE_PS_TP_DECL>
    template <class T>
    inline void base_pseries<__PIRANHA_BASE_PS_TP>::generic_builder(const T &x)
  {
    term_type term = term_type(cf_type(x),trig_type());
    if (!term.cf().is_insertable(cf_width()))
    {
      std::cout << "Warning: too many arguments in ctor from coefficient, building null series." << std::endl;
      return;
    }
    insert_check_positive(term);
  }
}

#endif
