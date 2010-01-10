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
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/tuple/tuple.hpp>
#include <utility>
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
		// Convert to max_fast_int.
		max_fast_int to_max_fast_int(const T &x) const
		{
			return boost::numeric_cast<max_fast_int>(x);
		}
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
		// Convert to max_fast_int.
		// NOTE: here we could benefit a lot from a nice implementation of numeric_cast and numeric traits...
		max_fast_int to_max_fast_int(const mp_rational &q) const
		{
			return boost::lexical_cast<max_fast_int>(q * m_lcd);
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
		// Decoding tuple: vectors of pairs of max_fast_ints, num and den resp. of the decoding formula.
		typedef boost::tuples::cons<std::vector<std::pair<max_fast_int,max_fast_int> >,
			typename cm_tuple_impl<typename Series::term_type::cf_type,N - 1>::type_decoding_tuple> type_decoding_tuple;
	};

	template <class Series>
	struct cm_tuple_impl<Series,0> {
		typedef boost::tuples::null_type type_minmax;
		typedef boost::tuples::null_type type_mp_minmax;
		typedef boost::tuples::null_type type_max_fast_int_minmax;
		typedef boost::tuples::null_type type_value_handler;
		typedef boost::tuples::null_type type_mp_coding_tuple;
		typedef boost::tuples::null_type type_coding_tuple;
		typedef boost::tuples::null_type type_decoding_tuple;
	};

	template <class Series>
	struct cm_tuple {
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_minmax type_minmax;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_mp_minmax type_mp_minmax;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_max_fast_int_minmax type_max_fast_int_minmax;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_value_handler type_value_handler;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_mp_coding_tuple type_mp_coding_tuple;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_coding_tuple type_coding_tuple;
		typedef typename cm_tuple_impl<Series,Series::echelon_level + 1>::type_decoding_tuple type_decoding_tuple;
	};

	template <class Series, class Tuple>
	struct cm_init_vector_tuples_impl {
		template <class ArgsTuple>
		static void run(const ArgsTuple &args_tuple, Tuple &t)
		{
			t.get_head().resize(args_tuple.template get<Series::term_type::key_type::position>().size());
			cm_init_vector_tuples_impl<typename Series::term_type::cf_type,typename Tuple::tail_type>::run(args_tuple,t.get_tail());
		}
	};

	template <class Series>
	struct cm_init_vector_tuples_impl<Series,boost::tuples::null_type>
	{
		template <class ArgsTuple>
		static void run(const ArgsTuple &, const boost::tuples::null_type &) {}
	};

	// Initialise the tuples-of-vectors types in coded multiplier with proper vector sizes.
	template <class Series, class Tuple, class ArgsTuple>
	inline void cm_init_vector_tuple(Tuple &t, const ArgsTuple &args_tuple)
	{
		p_static_check(boost::tuples::length<Tuple>::value == Series::echelon_level + 1,"");
		cm_init_vector_tuples_impl<Series,Tuple>::run(args_tuple,t);
	}

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
				tmp1.assign(boost::lexical_cast<mp_integer>(minmax1.get_head()[i].lower()),
					boost::lexical_cast<mp_integer>(minmax1.get_head()[i].upper()));
				tmp2.assign(boost::lexical_cast<mp_integer>(minmax2.get_head()[i].lower()),
					boost::lexical_cast<mp_integer>(minmax2.get_head()[i].upper()));
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
		// NOTE: maybe here it can be done better like below in cm_build_decoding_tuple.
		p_static_check(boost::tuples::length<MpCt>::value == boost::tuples::length<MpMinMaxTuple>::value,"");
		mp_integer const *prev_value = 0;
		boost::numeric::interval<mp_integer> const *prev_interval = 0;
		cm_mp_ct_impl<MpCt>::run(prev_value,prev_interval,mp_ct,mp_gt);
	}

	template <class Tuple1, class Tuple2>
	struct tuple_vector_dot_impl {
		p_static_check(boost::tuples::length<Tuple1>::value == boost::tuples::length<Tuple2>::value,"");
		typedef typename Tuple1::head_type::value_type value_type;
		typedef typename Tuple1::head_type::size_type size_type;
		p_static_check((boost::is_same<size_type,typename Tuple2::head_type::size_type>::value),"");
		static void run(const Tuple1 &t1, const Tuple2 &t2, value_type &retval)
		{
			piranha_assert(t1.get_head().size() == t2.get_head().size());
			for (size_type i = 0; i < t1.get_head().size(); ++i) {
				retval += t1.get_head()[i] * t2.get_head()[i];
			}
			tuple_vector_dot_impl<typename Tuple1::tail_type, typename Tuple2::tail_type>::run(
				t1.get_tail(),t2.get_tail(),retval);
		}
	};

	template <>
	struct tuple_vector_dot_impl<boost::tuples::null_type,boost::tuples::null_type> {
		template <class T>
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &, const T &) {}
	};

	// Calculate the dot product between two tuples of vectors. The result is assumed to be of the same type
	// of the values held in the head of the tuple.
	template <class Tuple1, class Tuple2>
	inline void tuple_vector_dot(const Tuple1 &t1, const Tuple2 &t2, typename Tuple1::head_type::value_type &retval)
	{
		tuple_vector_dot_impl<Tuple1,Tuple2>::run(t1,t2,retval);
	}

	template <class MpTuple>
	struct cm_mp_tuple_downcast_impl {
		template <class FastTuple>
		static void run(const MpTuple &mp_tuple, FastTuple &fast_tuple)
		{
			run_impl(mp_tuple.get_head(),fast_tuple.get_head());
			cm_mp_tuple_downcast_impl<typename MpTuple::tail_type>::run(mp_tuple.get_tail(),fast_tuple.get_tail());
		}
		static void run_impl(const std::vector<mp_integer> &mp_vector, std::vector<max_fast_int> &fast_vector)
		{
			typedef std::vector<mp_integer>::size_type size_type;
			p_static_check((boost::is_same<size_type,std::vector<max_fast_int>::size_type>::value),"");
			piranha_assert(mp_vector.size() == fast_vector.size());
			const size_type size = mp_vector.size();
			for (size_type i = 0; i < size; ++i) {
				fast_vector[i] = boost::lexical_cast<max_fast_int>(mp_vector[i]);
			}
		}
		static void run_impl(const std::vector<boost::numeric::interval<mp_integer> > &mp_vector,
			std::vector<boost::numeric::interval<max_fast_int> > &fast_vector)
		{
			typedef std::vector<boost::numeric::interval<mp_integer> >::size_type size_type;
			p_static_check((boost::is_same<size_type,std::vector<boost::numeric::interval<max_fast_int> >::size_type>::value),"");
			piranha_assert(mp_vector.size() == fast_vector.size());
			const size_type size = mp_vector.size();
			for (size_type i = 0; i < size; ++i) {
				fast_vector[i].assign(boost::lexical_cast<max_fast_int>(mp_vector[i].lower()),
				boost::lexical_cast<max_fast_int>(mp_vector[i].upper()));
			}
		}
	};

	template <>
	struct cm_mp_tuple_downcast_impl<boost::tuples::null_type> {
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

	// Downcast a tuple of multiprecision integer (interval) vectors to tuple of max_fast_int integer
	// (interval) vectors.
	template <class MpTuple, class FastTuple>
	inline void cm_mp_tuple_downcast(const MpTuple &mp_tuple, FastTuple &fast_tuple)
	{
		cm_mp_tuple_downcast_impl<MpTuple>::run(mp_tuple,fast_tuple);
	}

	// Check whether operations tuple contains at least one subtraction.
	template <class OpTuple>
	struct op_has_sub {
		static const bool value = (!OpTuple::head_type::value) || op_has_sub<typename OpTuple::tail_type>::value;
	};

	template <>
	struct op_has_sub<boost::tuples::null_type> {
		static const bool value = false;
	};

	template <class OpTuple>
	struct cm_code_impl2 {
		template <class CodingTuple, class Cf, class VhTuple>
		static void run(const CodingTuple &ct, const Cf &cf, const VhTuple &vh_tuple, max_fast_int &retval1, max_fast_int &retval2)
		{
			p_static_check((boost::is_same<typename CodingTuple::head_type::value_type,max_fast_int>::value),"");
			piranha_assert(cf.length() == 1);
			piranha_assert(cf.begin()->m_key.size() <= ct.get_head().size());
			typedef typename Cf::term_type::key_type::size_type size_type;
			max_fast_int tmp = 0;
			for (size_type i = 0; i < cf.begin()->m_key.size(); ++i) {
				tmp = ct.get_head()[i] * vh_tuple.get_head().to_max_fast_int(cf.begin()->m_key[i]);
				retval1 += tmp;
				if (!OpTuple::head_type::value) {
					retval2 -= tmp;
				} else {
					retval2 += tmp;
				}
			}
			cm_code_impl2<typename OpTuple::tail_type>::run(ct.get_tail(),cf.begin()->m_cf,vh_tuple.get_tail(),retval1,retval2);
		}
	};

	template <>
	struct cm_code_impl2<boost::tuples::null_type> {
		template <class Cf>
		static void run(const boost::tuples::null_type &, const Cf &, const boost::tuples::null_type &, const max_fast_int &, const max_fast_int &) {}
	};

	template <class OpTuple>
	struct cm_code_impl1 {
		template <class CodingTuple, class Term, class VhTuple>
		static void run(const CodingTuple &ct, const Term &term, const VhTuple &vh_tuple, max_fast_int &retval1, max_fast_int &retval2)
		{
			p_static_check((boost::is_same<typename CodingTuple::head_type::value_type,max_fast_int>::value),"");
			piranha_assert(term.m_key.size() <= ct.get_head().size());
			typedef typename Term::key_type::size_type size_type;
			max_fast_int tmp = 0;
			// NOTE: again the assumption that the sizes of vector and key are compatible. Need to sort this out...
			for (size_type i = 0; i < term.m_key.size(); ++i) {
				tmp = ct.get_head()[i] * vh_tuple.get_head().to_max_fast_int(term.m_key[i]);
				retval1 += tmp;
				if (!OpTuple::head_type::value) {
					retval2 -= tmp;
				} else {
					retval2 += tmp;
				}
			}
			cm_code_impl2<typename OpTuple::tail_type>::run(ct.get_tail(),term.m_cf,vh_tuple.get_tail(),retval1,retval2);
		}
	};

	// Code term using provided operations tuple type, coding tuple, value handler tuple and appending result to vector of codes v_codes.
	template <class OpTuple, class CodingTuple, class Term, class VhTuple>
	inline void cm_code(const CodingTuple &ct, const Term &term, const VhTuple &vh_tuple, max_fast_int &retval1, max_fast_int &retval2)
	{
		p_static_check(boost::tuples::length<OpTuple>::value == boost::tuples::length<CodingTuple>::value,"");
		p_static_check(boost::tuples::length<OpTuple>::value == boost::tuples::length<VhTuple>::value,"");
		// Initialise return values.
		retval1 = 0;
		retval2 = 0;
		cm_code_impl1<OpTuple>::run(ct,term,vh_tuple,retval1,retval2);
	}

	template <class DecodingTuple>
	struct cm_build_decoding_tuple_impl {
		template <class MinMaxTuple>
		static void run(max_fast_int *prev_range, DecodingTuple &dt, const MinMaxTuple &minmax)
		{
			typedef typename DecodingTuple::head_type::size_type size_type;
			const size_type size = boost::numeric_cast<size_type>(minmax.get_head().size());
			dt.get_head().resize(size);
			for (size_type i = 0; i < size; ++i) {
				dt.get_head()[i].first = (boost::numeric::width(minmax.get_head()[i]) + 1) * (*prev_range);
				dt.get_head()[i].second = *prev_range;
				prev_range = &dt.get_head()[i].first;
			}
			cm_build_decoding_tuple_impl<typename DecodingTuple::tail_type>::run(prev_range,dt.get_tail(),minmax.get_tail());
		}
	};

	template <>
	struct cm_build_decoding_tuple_impl<boost::tuples::null_type> {
		static void run(max_fast_int *, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

	// Build decoding tuple from global minmax representation.
	template <class DecodingTuple, class MinMaxTuple>
	inline void cm_build_decoding_tuple(DecodingTuple &dt, const MinMaxTuple &minmax)
	{
		p_static_check(boost::tuples::length<DecodingTuple>::value == boost::tuples::length<MinMaxTuple>::value,"");
		max_fast_int prev_range = 1;
		cm_build_decoding_tuple_impl<DecodingTuple>::run(&prev_range,dt,minmax);
	}

	template <class DecodingTuple>
	struct cm_decode_impl2 {
		template <class FinalCf, class MinMaxTuple, class Cf, class VhTuple, class ArgsTuple>
		static void run(const FinalCf &final_cf, const DecodingTuple &dt, const MinMaxTuple &gr, Cf &cf,
			const VhTuple &vh_tuple, const max_fast_int &code, const ArgsTuple &args_tuple)
		{
			piranha_assert(cf.empty());
			typedef typename Cf::term_type term_type;
			typedef typename term_type::key_type::size_type size_type;
			piranha_assert(dt.get_head().size() ==
				boost::numeric_cast<typename DecodingTuple::head_type::size_type>(args_tuple.template get<term_type::key_type::position>().size()));
			piranha_assert(dt.get_head().size() ==
				boost::numeric_cast<typename DecodingTuple::head_type::size_type>(gr.get_head().size()));
			const size_type size = boost::numeric_cast<size_type>(dt.get_head().size());
			// This term is going to be inserted into the coefficient series.
			term_type term;
			term.m_key.resize(size);
			for (size_type i = 0; i < size; ++i) {
				const max_fast_int tmp = (code % dt.get_head()[i].first) / dt.get_head()[i].second +
					gr.get_head()[i].lower();
				vh_tuple.get_head().post_decode(term.m_key[i],tmp);
			}
			// Next recursion will operate on the term-to-be-inserted's coefficient.
			cm_decode_impl2<typename DecodingTuple::tail_type>::run(final_cf,dt.get_tail(),gr.get_tail(),term.m_cf,vh_tuple.get_tail(),code,args_tuple);
			cf.insert(term,args_tuple);
		}
	};

	template <>
	struct cm_decode_impl2<boost::tuples::null_type> {
		template <class Cf, class ArgsTuple>
		static void run(const Cf &final_cf, const boost::tuples::null_type &, const boost::tuples::null_type &, Cf &cf,
			const boost::tuples::null_type &, const max_fast_int &, const ArgsTuple &)
		{
			// Last iteration simply assigns the final coefficient.
			cf = final_cf;
		}
	};

	template <class DecodingTuple>
	struct cm_decode_impl1 {
		// NOTE: here the code is not the shifted one, it is supposed to have been shifted in
		// the outside calling function.
		template <class FinalCf, class MinMaxTuple, class Term, class VhTuple, class ArgsTuple>
		static void run(const FinalCf &final_cf, const DecodingTuple &dt, const MinMaxTuple &gr, Term &term, const VhTuple &vh_tuple,
			const max_fast_int &code, const ArgsTuple &args_tuple)
		{
			piranha_assert(code >= 0);
			typedef typename Term::key_type::size_type size_type;
			// Make sure the sizes of the current tuples' elements are consistent.
			// (De)coding tuple sizes should always be the same as args tuple's.
			piranha_assert(dt.get_head().size() ==
				boost::numeric_cast<typename DecodingTuple::head_type::size_type>(args_tuple.template get<Term::key_type::position>().size()));
			piranha_assert(dt.get_head().size() ==
				boost::numeric_cast<typename DecodingTuple::head_type::size_type>(gr.get_head().size()));
			const size_type size = boost::numeric_cast<size_type>(dt.get_head().size());
			term.m_key.resize(size);
			for (size_type i = 0; i < size; ++i) {
				const max_fast_int tmp = (code % dt.get_head()[i].first) / dt.get_head()[i].second +
					gr.get_head()[i].lower();
				vh_tuple.get_head().post_decode(term.m_key[i],tmp);
			}
			cm_decode_impl2<typename DecodingTuple::tail_type>::run(final_cf,dt.get_tail(),gr.get_tail(),term.m_cf,
				vh_tuple.get_tail(),code,args_tuple);
		}
	};

	// Decode given code into term. final_cf is the last coefficient in the echelon hierarchy, from outwards to inwards.
	template <class FinalCf, class DecodingTuple, class MinMaxTuple, class Term, class VhTuple, class ArgsTuple>
	inline void cm_decode(const FinalCf &final_cf, const DecodingTuple &dt, const MinMaxTuple &gr, Term &term, const VhTuple &vh_tuple,
		const max_fast_int &code, const max_fast_int &min_code, const ArgsTuple &args_tuple)
	{
		p_static_check(boost::tuples::length<DecodingTuple>::value == boost::tuples::length<VhTuple>::value,"");
		p_static_check(boost::tuples::length<DecodingTuple>::value == boost::tuples::length<MinMaxTuple>::value,"");
		cm_decode_impl1<DecodingTuple>::run(final_cf,dt,gr,term,vh_tuple,code - min_code,args_tuple);
	}
}

#endif
