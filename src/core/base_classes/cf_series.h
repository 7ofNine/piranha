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

#ifndef PIRANHA_CF_SERIES_H
#define PIRANHA_CF_SERIES_H

#include <complex>
#include <iostream>
#include <string>

#include "base_series.h"
#include "../integer_typedefs.h" // For ctor macros
#include "../exceptions.h"
#include "../type_traits.h" // For specialisation of type trait

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_CF_SERIES_TP_DECL class Derived
#define __PIRANHA_CF_SERIES_TP Derived

namespace piranha
{
	/// Toolbox for using a series as a coefficient in another series.
	/**
	 * Intended to be inherited by piranha::base_series.
	 */
	template <__PIRANHA_CF_SERIES_TP_DECL>
	class cf_series
	{
		protected:
			// Proxy implementation to be refined in derived class.
			// Templated this way to cope with complex multiplication.
			template <class Series>
			class reference_proxy
			{
					friend class cf_series;
				public:
					reference_proxy(const Series &d):m_ptr(&d) {}
				protected:
					const Series	*m_ptr;
			};
		public:
			template <class ArgsTuple>
			void print_plain(std::ostream &, const ArgsTuple &) const;
			template <class ArgsTuple>
			bool is_insertable(const ArgsTuple &) const;
			template <class ArgsTuple>
			bool needs_padding(const ArgsTuple &) const;
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const;
			template <class ArgsTuple>
			void pad_right(const ArgsTuple &);
			template <class ArgsTuple>
			void invert_sign(const ArgsTuple &);
			void swap(Derived &);
			template <class ArgsTuple, class Layout>
			void apply_layout(const ArgsTuple &, const Layout &);
			template <class TrimFlags>
			void trim_test(TrimFlags &) const;
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &, const ArgsTuple &) const;
			// Interaction with proxy.
			template <class Series, class ArgsTuple>
			Derived &mult_by(reference_proxy<Series> , const ArgsTuple &);
		protected:
			template <class ArgsTuple>
			void construct_from_string(const std::string &, const ArgsTuple &);
	};

	// Useful macro for ctors in coefficient series.
	// TODO: maybe we can call these base_series ctors and use them in named_series ctors macro too?
#define CF_SERIES_CTORS(series_name) \
	series_name() {} \
	template <class ArgsTuple> \
	explicit series_name(const std::string &s, const ArgsTuple &args_tuple) \
	{ \
		cf_ancestor::construct_from_string(s,args_tuple); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const piranha::max_fast_int &n, const ArgsTuple &a) \
	{ \
		base_ancestor::construct_from_number(n,a); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const double &x, const ArgsTuple &a) \
	{ \
		base_ancestor::construct_from_number(x,a); \
	} \
	series_name(proxy p) { \
		*this=(*p.m_ptr); \
	} \
	series_name &operator=(proxy p) { \
		*this=(*p.m_ptr); \
		return *this; \
	}

#define COMPLEX_CF_SERIES_CTORS(complex_toolbox) \
	template <class ArgsTuple> \
	explicit complex(const complex<piranha::max_fast_int> &cn, const ArgsTuple &a) { \
		base_ancestor::construct_from_number(cn,a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const complex<double> &cx, const ArgsTuple &a) { \
		base_ancestor::construct_from_number(cx,a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const value_type &r, const ArgsTuple &a) { \
		complex_toolbox::construct_from_real(r,a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const value_type &r, const value_type &i, const ArgsTuple &a) { \
		complex_toolbox::construct_from_real_imag(r, i, a); \
	}

#define CF_SERIES_TERM(term_name,separator) term_name<Cf,Key,separator,Allocator>
#define CF_SERIES_BASE_ANCESTOR(term_name,series_name,term_separator,separator) \
	piranha::base_series<CF_SERIES_TERM(term_name,term_separator),separator, \
	Allocator,E0_SERIES(series_name) >

#define COMPLEX_CF_SERIES_TERM(term_name,separator) term_name<std::complex<Cf>,Key,separator,Allocator>
#define COMPLEX_CF_SERIES_BASE_ANCESTOR(term_name,series_name,term_separator,separator) piranha::base_series< \
	COMPLEX_CF_SERIES_TERM(term_name,term_separator),separator, Allocator,COMPLEX_E0_SERIES(series_name) >
}

#include "cf_series_io.h"
#include "cf_series_manip.h"
#include "cf_series_math.h"
#include "cf_series_probe.h"

#undef __PIRANHA_CF_SERIES_TP
#undef __PIRANHA_CF_SERIES_TP_DECL
#undef derived_const_cast
#undef derived_cast

#endif
