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

#ifndef PIRANHA_CODED_MULTIPLIER_MP_H
#define PIRANHA_CODED_MULTIPLIER_MP_H

#include <boost/lexical_cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"

// Meta-programming methods for coded multiplier.

namespace piranha {
	// Default value handler class.
	template <class T>
	struct cm_value_handler
	{
		// Assign to the minmax vector the values in the array key.
		template <class Key>
		void assign(std::vector<boost::numeric::interval<T> > &minmax, const Key &key)
		{
			p_static_check((boost::is_same<T,typename Key::value_type>::value),"");
			piranha_assert(key.size() <= minmax.size());
			for (typename Key::size_type i = 0; i < key.size(); ++i) {
				minmax[i] = key[i];
			}
		}
		// Test the values in array key and, if they sit outside the corresponding minmax intervals,
		// update the intervals to include them.
		template <class Key>
		void test(std::vector<boost::numeric::interval<T> > &minmax, const Key &key)
		{
			p_static_check((boost::is_same<T,typename Key::value_type>::value),"");
			piranha_assert(key.size() <= minmax.size());
			for (typename Key::size_type i = 0; i < key.size(); ++i) {
				if (key[i] < minmax[i].lower()) {
					minmax[i].assign(key[i],minmax[i].upper());
				} else if (key[i] > minmax[i].upper()) {
					minmax[i].assign(minmax[i].lower(),key[i]);
				}
			}
		}
		// Harmonization of two vh won't do anything by default.
		static void harmonize(cm_value_handler &, std::vector<boost::numeric::interval<T> > &,
			cm_value_handler &, std::vector<boost::numeric::interval<T> > &) {}
	};

	// Value handler class for rationals.
	// Similar to the default one, apart from the fact that all rationals will have denominator equal to one and
	// we keep trace of the global lcd.
	template <>
	struct cm_value_handler<mp_rational>
	{
		cm_value_handler(): m_lcd(1),m_old_lcd(1),m_tmp(0) {}
		template <class Key>
		void assign(std::vector<boost::numeric::interval<mp_rational> > &minmax, const Key &key)
		{
			p_static_check((boost::is_same<mp_rational,typename Key::value_type>::value),"");
			piranha_assert(key.size() <= minmax.size());
			for (typename Key::size_type i = 0; i < key.size(); ++i) {
				compute_new_lcd_and_update(minmax,key[i].get_den());
				// Assign current value.
				minmax[i] = key[i];
				// Assign tmp value.
				m_tmp = m_lcd;
				// Scale interval down to the lcd.
				minmax[i] *= m_tmp;
				piranha_assert(minmax[i].lower().get_den() == 1 && minmax[i].upper().get_den() == 1);
			}
		}
		template <class Key>
		void test(std::vector<boost::numeric::interval<mp_rational> > &minmax, const Key &key)
		{
			p_static_check((boost::is_same<mp_rational,typename Key::value_type>::value),"");
			piranha_assert(key.size() <= minmax.size());
			for (typename Key::size_type i = 0; i < key.size(); ++i) {
				// Compute the new lcd and update the current minmax_vector.
				compute_new_lcd_and_update(minmax,key[i].get_den());
				// Assign tmp value.
				m_tmp = key[i];
				// Scale interval down to the new lcd.
				m_tmp *= m_lcd;
				piranha_assert(m_tmp.get_den() == 1);
				if (m_tmp < minmax[i].lower()) {
					minmax[i].assign(m_tmp,minmax[i].upper());
				} else if (m_tmp > minmax[i].upper()) {
					minmax[i].assign(minmax[i].lower(),m_tmp);
				}
				piranha_assert(minmax[i].lower().get_den() == 1 && minmax[i].upper().get_den() == 1);
			}
		}
		// Harmonization here implies that we scale down the two vh to the same global lcd, and we update the corresponding
		// minmax vectors.
		static void harmonize(cm_value_handler &vh1, std::vector<boost::numeric::interval<mp_rational> > &minmax1,
			cm_value_handler &vh2, std::vector<boost::numeric::interval<mp_rational> > &minmax2)
		{
			vh1.compute_new_lcd_and_update(minmax1,vh2.m_lcd);
			vh2.compute_new_lcd_and_update(minmax2,vh1.m_lcd);
		}
		// Store current lcd into m_old_lcd, compute lcd between current lcd and input argument value,
		// store it into m_lcd and update the minmax vector to take into account the new lcd if needed.
		void compute_new_lcd_and_update(std::vector<boost::numeric::interval<mp_rational> > &minmax, const mp_integer &value)
		{
			// NOTE: here we assume that std::vector's size is greater than or equal to the corresponding key size.
			// This is probably "practically" true everywhere, but we should put a static assert on size_type
			// into key arrays, just in case. Or maybe wherever we make this assumption? Mmmh... Maybe we should derive
			// key arrays from common base class and put static checks in there, like range, unsignedness, etc.
			typedef std::vector<boost::numeric::interval<mp_rational> >::size_type size_type;
			// Store the old lcd.
			m_old_lcd = m_lcd;
			// Update the current lcd.
			m_lcd.lcm(m_lcd,value);
			// If the lcd changed, we need to update all the values
			// of the interval vector.
			if (m_old_lcd != m_lcd) {
				const size_type size = minmax.size();
				for (size_type i = 0; i < size; ++i) {
					m_tmp = m_lcd;
					minmax[i] *= m_tmp;
					m_tmp = m_old_lcd;
					minmax[i] /= m_tmp;
					piranha_assert(minmax[i].lower().get_den() == 1 && minmax[i].upper().get_den() == 1);
				}
			}
		}
		// Lowest common denominator.
		mp_integer	m_lcd;
		// Old value of the lcd.
		mp_integer	m_old_lcd;
		// Useful temporary rational value.
		mp_rational	m_tmp;
	};

	// Define a type to hold the min/max values of array keys in series.
	template <class Series, int N>
	struct cm_tuple_impl {
		p_static_check(N > 0,"");
		// minmax type, to be used for limits of input series.
		typedef typename Series::term_type::key_type::value_type value_type;
		typedef boost::tuples::cons<std::vector<boost::numeric::interval<value_type> >,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_minmax> type_minmax;
		// mp_integer minmax type, to be used when calculating global limits.
		typedef boost::tuples::cons<std::vector<boost::numeric::interval<mp_integer> >,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_mp_minmax> type_mp_minmax;
		// max_fast_int minmax type.
		typedef boost::tuples::cons<std::vector<boost::numeric::interval<max_fast_int> >,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_max_fast_int_minmax> type_max_fast_int_minmax;
		// value_handler tuple.
		typedef boost::tuples::cons<cm_value_handler<value_type>,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_value_handler> type_value_handler;
		// Multi-precision coding tuple.
		typedef boost::tuples::cons<std::vector<mp_integer>,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_mp_coding_tuple> type_mp_coding_tuple;
		// Coding tuple.
		typedef boost::tuples::cons<std::vector<max_fast_int>,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_coding_tuple> type_coding_tuple;
	};

	template <class Series>
	struct cm_tuple_impl<Series,0> {
		typedef boost::tuples::null_type type_minmax;
		typedef boost::tuples::null_type type_mp_minmax;
		typedef boost::tuples::null_type type_max_fast_int_minmax;
		typedef boost::tuples::null_type type_value_handler;
		typedef boost::tuples::null_type type_mp_coding_tuple;
		typedef boost::tuples::null_type type_coding_tuple;
	};

	template <class Series>
	struct cm_tuple {
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_minmax type_minmax;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_mp_minmax type_mp_minmax;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_max_fast_int_minmax type_max_fast_int_minmax;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_value_handler type_value_handler;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_mp_coding_tuple type_mp_coding_tuple;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_coding_tuple type_coding_tuple;
	};

	// Initialise the tuples-of-vectors types in coded multiplier with proper vector sizes.
	template <class Series, class Tuple>
	struct cm_init_vector_tuples {
		template <class ArgsTuple>
		static void run(const ArgsTuple &args_tuple, Tuple &t)
		{
			t.get_head().resize(args_tuple.template get<Series::term_type::key_type::position>().size());
			cm_init_vector_tuples<typename Series::term_type::cf_type,typename Tuple::tail_type>::run(args_tuple,t.get_tail());
		}
	};

	template <class Series>
	struct cm_init_vector_tuples<Series,boost::tuples::null_type>
	{
		template <class ArgsTuple>
		static void run(const ArgsTuple &, const boost::tuples::null_type &) {}
	};

	template <class MinMaxTuple>
	struct cm_minmax2 {
		template <class Series, class ValueHandlerTuple>
		static void run_init(const Series &s, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			p_static_check(boost::tuples::length<MinMaxTuple>::value == boost::tuples::length<ValueHandlerTuple>::value,"");
			piranha_assert(s.length() == 1);
			// NOTE: here key size could be less than the size of the vector in this tuple position,
			// hence the importance of having the vector default-initialised to zero. This happens because
			// we allow multiplication by series with fewer arguments.
			vh_tuple.get_head().assign(minmax_tuple.get_head(),s.begin()->m_key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_init(s.begin()->m_cf,minmax_tuple.get_tail(),vh_tuple.get_tail());
		}
		template <class Series, class ValueHandlerTuple>
		static void run_test(const Series &s, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			piranha_assert(s.length() == 1);
			vh_tuple.get_head().test(minmax_tuple.get_head(),s.begin()->m_key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_test(s.begin()->m_cf,minmax_tuple.get_tail(),vh_tuple.get_tail());
		}
	};

	template <>
	struct cm_minmax2<boost::tuples::null_type> {
		template <class Series>
		static void run_init(const Series &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
		template <class Series>
		static void run_test(const Series &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	template <class MinMaxTuple>
	struct cm_minmax {
		// Initialise minmax tuples with a term.
		template <class Term, class ValueHandlerTuple>
		static void run_init(const Term &term, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			p_static_check(boost::tuples::length<MinMaxTuple>::value == boost::tuples::length<ValueHandlerTuple>::value,"");
			vh_tuple.get_head().assign(minmax_tuple.get_head(),term.m_key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_init(term.m_cf,minmax_tuple.get_tail(),vh_tuple.get_tail());
		}
		// Test minmax tuple with term's array representation, and, if needed, update min and max to include it.
		template <class Term, class ValueHandlerTuple>
		static void run_test(const Term &term, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			vh_tuple.get_head().test(minmax_tuple.get_head(),term.m_key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_test(term.m_cf,minmax_tuple.get_tail(),vh_tuple.get_tail());
		}
	};

	// Calculate the minmax values of the global representation, using multiprecision integers.
	template <class OpTuple>
	struct cm_global_minmax
	{
		template <class MinMaxTuple, class ValueHandlerTuple, class MpMinMaxTuple>
		static void run(
			MinMaxTuple &minmax1, ValueHandlerTuple &vh1,
			MinMaxTuple &minmax2, ValueHandlerTuple &vh2,
			MpMinMaxTuple &global_minmax)
		{
			// Make sure sizes are consistent.
			p_static_check(boost::tuples::length<OpTuple>::value == boost::tuples::length<MinMaxTuple>::value,"");
			p_static_check(boost::tuples::length<OpTuple>::value == boost::tuples::length<MpMinMaxTuple>::value,"");
			p_static_check(boost::tuples::length<OpTuple>::value == boost::tuples::length<ValueHandlerTuple>::value,"");
			p_static_check((boost::is_same<typename MinMaxTuple::head_type::size_type,
				typename MpMinMaxTuple::head_type::size_type>::value),"");
			piranha_assert(minmax1.get_head().size() == minmax2.get_head().size() && 
				global_minmax.get_head().size() == minmax2.get_head().size());
			// Harmonize the value handler tuples.
			ValueHandlerTuple::head_type::harmonize(vh1.get_head(),minmax1.get_head(),vh2.get_head(),minmax2.get_head());
			boost::numeric::interval<mp_integer> tmp1, tmp2;
			for (typename MinMaxTuple::head_type::size_type i = 0; i < minmax1.get_head().size(); ++i) {
				// Transform the current intervals into mp intervals.
				tmp1.assign(mp_integer(boost::lexical_cast<std::string>(minmax1.get_head()[i].lower())),
					mp_integer(boost::lexical_cast<std::string>(minmax1.get_head()[i].upper())));
				tmp2.assign(mp_integer(boost::lexical_cast<std::string>(minmax2.get_head()[i].lower())),
					mp_integer(boost::lexical_cast<std::string>(minmax2.get_head()[i].upper())));
				// Addition.
				global_minmax.get_head()[i] = boost::numeric::hull(
					boost::numeric::hull(
						tmp1,tmp2
					),
						tmp1 + tmp2
				);
				// Perform also subtraction, if requested.
				if (!OpTuple::head_type::value) {
					global_minmax.get_head()[i] = boost::numeric::hull(
						global_minmax.get_head()[i],
						tmp1 - tmp2
					);
				}
			}
			cm_global_minmax<typename OpTuple::tail_type>::run(minmax1.get_tail(),vh1.get_tail(),
				minmax2.get_tail(),vh2.get_tail(),global_minmax.get_tail());
		}
	};

	template <>
	struct cm_global_minmax<boost::tuples::null_type>
	{
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &,
			const boost::tuples::null_type &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

	template <class MpCt>
	struct cm_mp_ct_impl {
		template <class MpMinMaxTuple>
		static void run(mp_integer const *prev_value, boost::numeric::interval<mp_integer> const *prev_interval,
			MpCt &mp_ct, const MpMinMaxTuple &mp_gt)
		{
			p_static_check(boost::tuples::length<MpCt>::value == boost::tuples::length<MpMinMaxTuple>::value,"");
			// Assume same sizes and size types.
			piranha_assert(mp_ct.get_head().size() == mp_gt.get_head().size());
			typedef typename MpCt::head_type::size_type size_type;
			p_static_check((boost::is_same<size_type,typename MpMinMaxTuple::head_type::size_type>::value),"");
			for (size_type i = 0; i < mp_ct.get_head().size(); ++i) {
				if (prev_interval) {
					piranha_assert(prev_value);
					// Assign to the current value the previous one.
					mp_ct.get_head()[i] = *prev_value;
					// Multiply it by previous interval's width + 1.
					mp_integer width(boost::numeric::width(*prev_interval));
					width += 1;
					mp_ct.get_head()[i] *= width;
				} else {
					// This corresponds to the first calculated element of the ct.
					// It *must* be at the beginning of a vector.
					piranha_assert(i == 0);
					piranha_assert(!prev_value);
					mp_ct.get_head()[i] = 1;
				}
				// Assign previous value and interval.
				prev_interval = &mp_gt.get_head()[i];
				prev_value = &mp_ct.get_head()[i];
			}
			cm_mp_ct_impl<typename MpCt::tail_type>::run(prev_value,prev_interval,mp_ct.get_tail(),mp_gt.get_tail());
		}
	};

	template <>
	struct cm_mp_ct_impl<boost::tuples::null_type> {
		static void run(mp_integer const *prev_value, boost::numeric::interval<mp_integer> const *prev_interval,
			const boost::tuples::null_type &, const boost::tuples::null_type &)
		{
			// Make sure that when we reach here we have computed at least one component of the coding tuple.
			// Otherwise it means that we are multiplying series with zero arguments and we should not be here.
			piranha_assert(prev_value);
			piranha_assert(prev_interval);
			(void)prev_value;
			(void)prev_interval;
		}
	};

	// Compute coding tuple from tuple of multi-precision minmax vectors.
	template <class MpCt, class MpMinMaxTuple>
	void compute_mp_coding_tuple(MpCt &mp_ct, const MpMinMaxTuple &mp_gt)
	{
		p_static_check(boost::tuples::length<MpCt>::value == boost::tuples::length<MpMinMaxTuple>::value,"");
		mp_integer const *prev_value = 0;
		boost::numeric::interval<mp_integer> const *prev_interval = 0;
		cm_mp_ct_impl<MpCt>::run(prev_value,prev_interval,mp_ct,mp_gt);
	}

	// Calculate the dot product between two tuples of vectors. The result is assumed to be of the same type
	// of the values held in the head of the tuple.
	template <class Tuple1, class Tuple2>
	struct tuple_vector_dot {
		p_static_check(boost::tuples::length<Tuple1>::value == boost::tuples::length<Tuple2>::value,"");
		typedef typename Tuple1::head_type::value_type value_type;
		typedef typename Tuple1::head_type::size_type size_type;
		static void run(const Tuple1 &t1, const Tuple2 &t2, value_type &retval)
		{
			piranha_assert(t1.get_head().size() == t2.get_head().size());
			for (size_type i = 0; i < t1.get_head().size(); ++i) {
				retval += t1.get_head()[i] * t2.get_head()[i];
			}
			tuple_vector_dot<typename Tuple1::tail_type, typename Tuple2::tail_type>::run(
				t1.get_tail(),t2.get_tail(),retval);
		}
	};

	template <>
	struct tuple_vector_dot<boost::tuples::null_type,boost::tuples::null_type> {
		template <class T>
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &, const T &) {}
	};

	// To test whether a representation is viable or not, we need to test for the following things:
	// - mp_h_minmax must be in the max_fast_int range;
	// - mp_h_minmax's width must be in the max_fast_int range.
}

#endif