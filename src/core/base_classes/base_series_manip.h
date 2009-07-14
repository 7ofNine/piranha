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

#ifndef PIRANHA_BASE_SERIES_MANIP_H
#define PIRANHA_BASE_SERIES_MANIP_H

#include <boost/tuple/tuple.hpp>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "base_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Transform term.
	/**
	 * Non-specialized version, it will create a copy of the converted term.
	 */
	template <class Term2, class Term1>
	class term_converter
	{
		public:
			/// Constructor.
			template <class ArgsTuple>
			explicit term_converter(const Term2 &c, const ArgsTuple &a): result(c,a) {}
			/// Copy of the converted term.
			const Term1 result;
	};

	/// Specialized term converter.
	/**
	 * It will be invoked when the type to convert from is the same as the converted type. A reference
	 * to the convertee is stored inside the class.
	 */
	template <class Term2>
	class term_converter<Term2, Term2>
	{
		public:
			/// Constructor.
			template <class ArgsTuple>
			explicit term_converter(const Term2 &c, const ArgsTuple &): result(c) {}
			/// Reference to the converted term.
			const Term2 &result;
	};

	// TODO: update doc here.
	/// High-level insertion function.
	/**
	 * This function is used to insert terms into a series. It requires that the number of arguments
	 * of each element of the term is smaller or equal to the series',
	 * otherwise an assertion fails and the program aborts. base_pseries::merge_args,
	 * base_pseries::append_cf_args, base_pseries::append_trig_args, etc. can be used to add the needed arguments
	 * to the series.
	 *
	 * This function performs some checks and then calls ll_insert.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool CanonicalCheck, bool Sign, class Term2, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::insert(const Term2 &term_, const ArgsTuple &args_tuple)
	{
		term_converter<Term2, term_type> converted_term(term_, args_tuple);
		// Make sure the appropriate routines for the management of arguments have been called.
		piranha_assert(converted_term.result.m_cf.is_insertable(args_tuple) &&
			converted_term.result.m_key.is_insertable(args_tuple));
		term_type *new_term(0);
		if (unlikely(converted_term.result.m_cf.needs_padding(args_tuple) ||
			converted_term.result.m_key.needs_padding(args_tuple))) {
			new_term = term_type::allocator.allocate(1);
			term_type::allocator.construct(new_term, converted_term.result);
			new_term->m_cf.pad_right(args_tuple);
			new_term->m_key.pad_right(args_tuple);
		}
		if (CanonicalCheck) {
			if (!converted_term.result.is_canonical(args_tuple)) {
				if (new_term == 0) {
					new_term = term_type::allocator.allocate(1);
					term_type::allocator.construct(new_term, converted_term.result);
				}
				new_term->canonicalise(args_tuple);
			}
		}
		const term_type *insert_term(0);
		if (new_term) {
			insert_term = new_term;
		} else {
			insert_term = &converted_term.result;
		}
		ll_insert<Sign>(*insert_term, args_tuple);
		if (new_term) {
			term_type::allocator.destroy(new_term);
			term_type::allocator.deallocate(new_term, 1);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Term2, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::insert(const Term2 &term, const ArgsTuple &args_tuple)
	{
		insert<true, true>(term, args_tuple);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::insert_range(const Iterator &begin,
		const Iterator &end, const ArgsTuple &args_tuple) {
		for (Iterator it = begin; it != end; ++it) {
			insert(*it,args_tuple);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::const_iterator
	toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::find_term(const term_type &t) const
	{
		return m_container.find(t);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::ll_insert(const term_type &term, const ArgsTuple &args_tuple)
	{
		// TODO: think about moving this check higher in the stack of functions, we probably don't want to reach
		// _this_ point before checking for ignorability.
		if (term.m_cf.is_ignorable(args_tuple) || term.m_key.is_ignorable(args_tuple)) {
			return;
		}
		piranha_assert(term.m_cf.is_insertable(args_tuple) && term.m_key.is_insertable(args_tuple) &&
			!term.m_cf.needs_padding(args_tuple) && !term.m_key.needs_padding(args_tuple) && term.is_canonical(args_tuple));
		const_iterator it(find_term(term));
		if (it == end()) {
			// The term is NOT a duplicate, insert in the set.
			term_insert_new<Sign>(term, args_tuple);
		} else {
			// The term is in the set, hence an existing term will be modified.
			// Add or subtract according to request.
			if (Sign) {
				it->m_cf.add(term.m_cf, args_tuple);
			} else {
				it->m_cf.subtract(term.m_cf, args_tuple);
			}
			// Check if the new coefficient can be ignored.
			if (it->m_cf.is_ignorable(args_tuple)) {
				term_erase(it, args_tuple);
			}
		}
	}

	// Insert a new term into the series
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::term_insert_new(const term_type &term,
			const ArgsTuple &args_tuple)
	{
		std::pair<const_iterator, bool> res(m_container.insert(term));
		piranha_assert(res.second);
		if (!Sign) {
			res.first->m_cf.invert_sign(args_tuple);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::term_erase(const const_iterator &it,
			const ArgsTuple &)
	{
		m_container.erase(it);
	}

	/// Swap the terms with another series.
	/**
	 * All terms get swapped.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_swap(Derived &ps2)
	{
		piranha_assert(derived_cast != &ps2);
		m_container.swap(ps2.m_container);
	}

	/// Apply an arguments layout to all terms and insert them into retval.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Layout, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::apply_layout_to_terms(
		const Layout &l, Derived &retval, const ArgsTuple &args_tuple) const
	{
		const const_iterator it_f = end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			term_type term(*it);
			term.m_cf.apply_layout(l,args_tuple);
			term.m_key.apply_layout(l,args_tuple);
			retval.insert(term, args_tuple);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::trim_test_terms(TrimFlags &tf) const
	{
		const const_iterator it_f = end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			it->m_cf.trim_test(tf);
			it->m_key.trim_test(tf);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::trim_terms(const TrimFlags &tf, Derived &retval,
			const ArgsTuple &args_tuple) const
	{
		const const_iterator it_f = end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			retval.insert(
				term_type(typename term_type::cf_type(it->m_cf.trim(tf, args_tuple)),
					typename term_type::key_type(it->m_key.trim(tf, args_tuple))),
				args_tuple
			);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
	inline RetSeries toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_sub(const PosTuple &pos_tuple,
			SubCaches &sub_caches, const ArgsTuple &args_tuple) const
	{
		p_static_check((boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value),
			"Positional and arguments' tuples' lengths do not match.");
		RetSeries retval;
		const const_iterator it_f = end();
		for (const_iterator it = begin(); it != it_f; ++it) {
			RetSeries tmp = SubFunctor::template run<RetSeries>(it->m_cf,pos_tuple,sub_caches,args_tuple);
			// NOTICE: series multadd here?
			tmp.base_mult_by(SubFunctor::template run<RetSeries>(it->m_key,pos_tuple,sub_caches,args_tuple),args_tuple);
			retval.base_add(tmp,args_tuple);
		}
		return retval;
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::clear_terms()
	{
		m_container.clear();
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Series, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::base_split(std::vector<std::vector<Series> > &retval, const int &n, const ArgsTuple &args_tuple) const
	{
		piranha_assert(retval.empty());
		piranha_assert(n >= 0 && n < boost::tuples::length<ArgsTuple>::value);
		if (n == 0) {
			try {
				const std::vector<typename Derived::term_type const *> s(derived_const_cast->template get_sorted_series<Derived>(args_tuple));
				generic_base_split(retval,&(*s.begin()),&(*s.end()),args_tuple);
			} catch (const value_error &) {
				generic_base_split(retval,begin(),end(),args_tuple);
			}
		} else {
			if (!is_single_cf()) {
				piranha_throw(value_error,"cannot split up to the specified level: series is non-degenerate");
			}
			begin()->m_cf.split(retval,n - 1,args_tuple);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class Series, class ArgsTuple>
	inline void toolbox<base_series<__PIRANHA_BASE_SERIES_TP> >::generic_base_split(std::vector<std::vector<Series> > &retval, const Iterator &start,
		const Iterator &end, const ArgsTuple &args_tuple) const
	{
		for (Iterator it = start; it != end; ++it) {
			Series tmp_cf(Series::base_series_from_cf(it_getter<Iterator>::get(it)->m_cf,args_tuple));
			Series tmp_key(Series::base_series_from_key(it_getter<Iterator>::get(it)->m_key,args_tuple));
			std::vector<Series> tmp;
			tmp.reserve(2);
			tmp.push_back(tmp_cf);
			tmp.push_back(tmp_key);
			retval.push_back(tmp);
		}
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
