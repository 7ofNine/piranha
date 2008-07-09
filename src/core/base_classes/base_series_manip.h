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

#include <boost/type_traits/is_same.hpp> // For iterator type detection.
#include <utility>

#include "../config.h"

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
			explicit term_converter(const Term2 &c, const ArgsTuple &a): result(c, a) {}
			/// Copy of the converted term.
			const Term1	result;
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
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term_, const ArgsTuple &args_tuple)
	{
		term_converter<Term2, term_type> converted_term(term_, args_tuple);
		// Make sure the appropriate routines for the management of arguments have been called.
		p_assert(converted_term.result.is_insertable(args_tuple));
		term_type *new_term(0);
		if (unlikely(converted_term.result.needs_padding(args_tuple))) {
			new_term = term_type::allocator.allocate(1);
			term_type::allocator.construct(new_term, converted_term.result);
			new_term->pad_right(args_tuple);
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
		if (new_term == 0) {
			insert_term = &converted_term.result;
		} else {
			insert_term = new_term;
		}
		ll_insert<Sign>(*insert_term, args_tuple);
		if (new_term != 0) {
			term_type::allocator.destroy(new_term);
			term_type::allocator.deallocate(new_term, 1);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Term2, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term, const ArgsTuple &args_tuple)
	{
		insert<true, true>(term, args_tuple);
	}

	// This cannot be const because we use this in insertion function, hence we need a non const iterator.
	/// Find term.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PinpointIterator>
	inline PinpointIterator base_series<__PIRANHA_BASE_SERIES_TP>::find_term(const term_type &t)
	{
		p_static_check((boost::is_same<PinpointIterator, typename Derived::pinpoint_iterator>::value), "");
		return derived_cast->template nth_index<1>().find(t);
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::ll_insert(const term_type &term, const ArgsTuple &args_tuple)
	{
		typedef typename Derived::pinpoint_iterator pinpoint_iterator;
		// TODO: think about moving this check higher in the stack of functions, we probably don't want to reach
		// _this_ point before checking for ignorability.
		if (term.is_ignorable(args_tuple)) {
			return;
		}
		p_assert(term.is_insertable(args_tuple) && !term.needs_padding(args_tuple) && term.is_canonical(args_tuple));
		pinpoint_iterator it(find_term<pinpoint_iterator>(term));
		if (it == derived_const_cast->template nth_index<1>().end()) {
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
				term_erase<1>(it, args_tuple);
			}
		}
	}

	// Insert a new term into the series
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::term_insert_new(const term_type &term,
			const ArgsTuple &args_tuple)
	{
		typedef typename Derived::const_sorted_iterator const_sorted_iterator;
		std::pair<const_sorted_iterator, bool> res(derived_cast->template nth_index<0>().push_back(term));
		p_assert(res.second);
		if (!Sign) {
			res.first->m_cf.invert_sign(args_tuple);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <int N, class Iterator, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::term_erase(Iterator it,
			const ArgsTuple &)
	{
		derived_cast->template nth_index<N>().erase(it);
	}

	/// Swap the terms with another series.
	/**
	 * All terms get swapped.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::swap_terms(Derived &ps2)
	{
		p_assert(derived_cast != &ps2);
		derived_cast->m_container.swap(ps2.m_container);
	}

	/// Apply an arguments layout to all terms and insert them into retval.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple, class Layout>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::apply_layout_to_terms(
		const ArgsTuple &args_tuple, const Layout &l, Derived &retval) const
	{
		typedef typename Derived::const_sorted_iterator const_sorted_iterator;
		typedef typename Derived::term_type term_type;
		const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
		for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin();
				it != it_f; ++it) {
			term_type term(*it);
			term.m_cf.apply_layout(args_tuple, l);
			term.m_key.apply_layout(args_tuple, l);
			retval.insert(term, args_tuple);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::trim_test_terms(TrimFlags &tf) const
	{
		typedef typename Derived::const_sorted_iterator const_sorted_iterator;
		const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
		for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it) {
			it->m_cf.trim_test(tf);
			it->m_key.trim_test(tf);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline void base_series<__PIRANHA_BASE_SERIES_TP>::trim_terms(const TrimFlags &tf, Derived &retval,
			const ArgsTuple &args_tuple) const
	{
		typedef typename Derived::const_sorted_iterator const_sorted_iterator;
		typedef typename Derived::term_type term_type;
		typedef typename term_type::cf_type cf_type;
		typedef typename term_type::key_type key_type;
		const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
		for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it) {
			retval.insert(
				term_type(cf_type(it->m_cf.trim(tf, args_tuple)), key_type(it->m_key.trim(tf, args_tuple))),
				args_tuple
			);
		}
	}

	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class PosTuple, class SubSeries, class ArgsTuple>
	inline Derived base_series<__PIRANHA_BASE_SERIES_TP>::sub(const PosTuple &pos_tuple, const SubSeries &s,
		const ArgsTuple &args_tuple) const {
		p_static_check((boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value),
			"Positional and arguments' tuples' lengths do not match in base_series::sub.");
		typedef typename Derived::const_sorted_iterator const_sorted_iterator;
		Derived retval;
		const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
		for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it) {
			Derived tmp(it->m_cf.sub<Derived>(pos_tuple,s,args_tuple));
			tmp.mult_by(it->m_key.sub<Derived>(pos_tuple,s,args_tuple),args_tuple);
			retval.add(tmp,args_tuple);
		}
		return retval;
	}
}

#endif
