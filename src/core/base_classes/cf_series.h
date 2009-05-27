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

#include <boost/tuple/tuple.hpp> // For sub cache.
#include <complex>
#include <iostream>
#include <string>

#include "../mp.h"
#include "../type_traits.h"
#include "base_series.h"
#include "toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_CF_SERIES_TP_DECL class Term, class Derived
#define __PIRANHA_CF_SERIES_TP Term,Derived

namespace piranha
{
	template <__PIRANHA_CF_SERIES_TP_DECL>
	struct cf_series {};

	/// Toolbox for using a series as a coefficient in another series.
	/**
	 * Intended to be inherited by piranha::base_series.
	 */
	template <__PIRANHA_CF_SERIES_TP_DECL>
	class toolbox<cf_series<__PIRANHA_CF_SERIES_TP> >
	{
		public:
			typedef typename term_eval_type_determiner<Term>::type eval_type;
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct sub_cache_selector {
				typedef typename Derived::term_type::cf_type::
					template sub_cache_selector<SubSeries,typename Derived::term_type::key_type::
					template sub_cache_selector<SubSeries,SubCachesCons,ArgsTuple>::type,ArgsTuple>::type type;
			};
			template <class ArgsTuple>
			void print_plain(std::ostream &, const ArgsTuple &) const;
			template <class ArgsTuple>
			void print_pretty(std::ostream &, const ArgsTuple &) const;
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
			template <class Layout, class ArgsTuple>
			void apply_layout(const Layout &, const ArgsTuple &);
			template <class TrimFlags>
			void trim_test(TrimFlags &) const;
			template <class TrimFlags, class ArgsTuple>
			Derived trim(const TrimFlags &, const ArgsTuple &) const;
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &, SubCaches &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived &add(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &add(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &add(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &add(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &subtract(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &subtract(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &subtract(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &subtract(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &mult_by(const Derived &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &mult_by(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &mult_by(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &mult_by(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &divide_by(const double &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &divide_by(const mp_rational &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived &divide_by(const mp_integer &, const ArgsTuple &);
			template <class ArgsTuple>
			Derived pow(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			Derived pow(const mp_rational &, const ArgsTuple &) const;
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &, const ArgsTuple &) const;
			template <class ArgsTuple>
			eval_type eval(const double &, const ArgsTuple &) const;
			template <class ArgsTuple>
			double norm(const ArgsTuple &) const;
			bool operator==(const Derived &) const;
			bool operator!=(const Derived &) const;
			bool operator==(const double &) const;
			bool operator!=(const double &) const;
			bool operator==(const mp_rational &) const;
			bool operator!=(const mp_rational &) const;
			bool operator==(const mp_integer &) const;
			bool operator!=(const mp_integer &) const;
		protected:
			template <class ArgsTuple>
			void construct_from_string(const std::string &, const ArgsTuple &);
	};

	// Useful macro for ctors in coefficient series.
	// TODO: maybe we can call these base_series ctors and use them in named_series ctors macro too?
#define CF_SERIES_CTORS(series_name) \
	explicit series_name() {} \
	template <class ArgsTuple> \
	explicit series_name(const std::string &s, const ArgsTuple &args_tuple) \
	{ \
		this->construct_from_string(s,args_tuple); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const double &x, const ArgsTuple &a) \
	{ \
		this->construct_from_number(x,a); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const piranha::mp_rational &q, const ArgsTuple &a) \
	{ \
		this->construct_from_number(q,a); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const piranha::mp_integer &z, const ArgsTuple &a) \
	{ \
		this->construct_from_number(z,a); \
	}

#define COMPLEX_CF_SERIES_CTORS(real_series) \
	template <class ArgsTuple> \
	explicit complex(const complex<double> &cx, const ArgsTuple &a) { \
		this->construct_from_number(cx,a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const real_series &r, const ArgsTuple &a) { \
		this->base_construct_from_real(r,a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const real_series &r, const real_series &i, const ArgsTuple &a) { \
		this->base_construct_from_real_imag(r, i, a); \
	}

#define CF_SERIES_TERM(term_name,separator) term_name<Cf,Key,separator,Allocator>
#define CF_SERIES_BASE_ANCESTOR(term_name,series_name,term_separator,separator) \
	piranha::toolbox<piranha::base_series<CF_SERIES_TERM(term_name,term_separator),separator, \
	Allocator,E0_SERIES(series_name) > >

#define COMPLEX_CF_SERIES_TERM(term_name,separator) term_name<std::complex<Cf>,Key,separator,Allocator>
#define COMPLEX_CF_SERIES_BASE_ANCESTOR(term_name,series_name,term_separator,separator) piranha::toolbox<piranha::base_series< \
	COMPLEX_CF_SERIES_TERM(term_name,term_separator),separator, Allocator,COMPLEX_E0_SERIES(series_name) > >
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
