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
#include <vector>

#include "../mp.h"
#include "../type_traits.h"
#include "base_series.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)
#define __PIRANHA_CF_SERIES_TP_DECL class Term, class Derived
#define __PIRANHA_CF_SERIES_TP Term,Derived

// TODO: split this off like in base/name series.

namespace piranha
{
	/// Toolbox for using a series as a coefficient in another series.
	// what is the difference to NamedSeries?
	/**
	 * Intended to be inherited by piranha::BaseSeries.
	 */
	template <__PIRANHA_CF_SERIES_TP_DECL>
	class cf_series
	{
		public:

			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct SubstitutionCacheSelector {
				typedef typename Derived::TermType::CfType::
					template SubstitutionCacheSelector<SubSeries, typename Derived::TermType::KeyType::
					template SubstitutionCacheSelector<SubSeries, SubCachesCons, ArgsTuple>::type, ArgsTuple>::type type;
			};

			template <class ArgsTuple>
			void print_plain(std::ostream &, const ArgsTuple &) const;

			template <class ArgsTuple>
			void print_pretty(std::ostream &, const ArgsTuple &) const;

			template <class ArgsTuple>
			void print_tex(std::ostream &, const ArgsTuple &) const;

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

			template <class T, class ArgsTuple>
			Derived &add(const T &, const ArgsTuple &);

			template <class T, class ArgsTuple>
			Derived &subtract(const T &, const ArgsTuple &);

			template <class T, class ArgsTuple>
			Derived &mult_by(const T &, const ArgsTuple &);

			template <class T, class ArgsTuple>
			Derived &divideBy(const T &, const ArgsTuple &);

			template <class ArgsTuple>
			Derived pow(const double &, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived pow(const mp_rational &, const ArgsTuple &) const;

			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &, const ArgsTuple &) const;

			template <class ArgsTuple>
			typename TermEvalTypeDeterminer<Term>::type eval(const double &, const ArgsTuple &) const;

			template <class ArgsTuple>
			double norm(const ArgsTuple &) const;

			template <class T>
			bool operator==(const T &) const;

			template <class T>
			bool operator!=(const T &) const;

			template <class Series, class ArgsTuple>
			void split(std::vector<std::vector<Series> > &, const int &, const ArgsTuple &) const;

		protected:

			template <class ArgsTuple>
			void construct_from_string(const std::string &, const ArgsTuple &);
	};


	// Useful macro for ctors in coefficient series.
	// TODO: maybe we can call these BaseSeries ctors and use them in NamedSeries ctors macro too?
#define CF_SERIES_CTORS(series_name) \
	explicit series_name() {} \
	template <class ArgsTuple> \
	explicit series_name(const std::string &s, const ArgsTuple &argsTuple) \
	{ \
		this->construct_from_string(s, argsTuple); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const double &x, const ArgsTuple &a) \
	{ \
		*this = baseSeriesFromNumber(x, a); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const piranha::mp_rational &q, const ArgsTuple &a) \
	{ \
		*this = baseSeriesFromNumber(q, a); \
	} \
	template <class ArgsTuple> \
	explicit series_name(const piranha::mp_integer &z, const ArgsTuple &a) \
	{ \
		*this = baseSeriesFromNumber(z, a); \
	}


#define COMPLEX_CF_SERIES_CTORS(real_series) \
	template <class ArgsTuple> \
	explicit complex(const complex<double> &cx, const ArgsTuple &a) { \
		*this = baseSeriesFromNumber(cx, a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const real_series &r, const ArgsTuple &a) { \
		this->baseConstructFromReal(r, a); \
	} \
	template <class ArgsTuple> \
	explicit complex(const real_series &r, const real_series &i, const ArgsTuple &a) { \
		this->baseConstructFromRealImag(r, i, a); \
	}


#define CF_SERIES_TERM(term_name, separator) term_name<Cf, Key, separator, Allocator>
#define CF_SERIES_BASE_ANCESTOR(term_name, series_name, term_separator, separator) \
	piranha::BaseSeries<CF_SERIES_TERM(term_name, term_separator), separator, \
	Allocator,E0_SERIES(series_name)>

#define COMPLEX_CF_SERIES_TERM(term_name, separator) term_name<std::complex<Cf>, Key, separator, Allocator>
#define COMPLEX_CF_SERIES_BASE_ANCESTOR(term_name, series_name, term_separator, separator) piranha::BaseSeries< \
	COMPLEX_CF_SERIES_TERM(term_name, term_separator), separator, Allocator, COMPLEX_E0_SERIES(series_name)>
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
