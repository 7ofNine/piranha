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
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/tuple/tuple.hpp>
#include <utility>
#include <vector>

#include "../base_classes/base_series_tag.h"
#include "../config.h"
#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../mp.h"

// Meta-programming methods for coded multiplier.

namespace piranha {
	// Getter method for final coefficient of flattened series term.
	template <class Cf, class Enable = void>
	struct final_cf_getter_impl
	{
		static const Cf &run(const Cf &cf)
		{
			return cf;
		}
	};

	template <class CfSeries>
	struct final_cf_getter_impl<CfSeries,typename boost::enable_if<boost::is_base_of<BaseSeriesTag,CfSeries> >::type>
	{
		static const typename final_cf<CfSeries>::type &run(const CfSeries &cf_series)
		{
			PIRANHA_ASSERT(cf_series.length() == 1);
			return final_cf_getter_impl<typename CfSeries::TermType::CfType>::run(cf_series.begin()->cf);
		}
	};


	template <class Series>
	struct final_cf_getter
	{
		const typename final_cf<Series>::type &operator()(const typename Series::TermType *term) const
		{
			return final_cf_getter_impl<typename Series::TermType::CfType>::run(term->cf);
		}
	};


	// Generalised reverse lexicographic ordering implementation.
	template <class Term, class Enable = void>
	struct key_revlex_comparison_impl
	{
		static bool run(const Term *t1, const Term *t2)
		{
			return t1->key.revlex_comparison(t2->key);
		}
	};


	template <class Term>
	struct key_revlex_comparison_impl<Term,typename boost::enable_if<boost::is_base_of<BaseSeriesTag,typename Term::CfType> >::type>
	{
		static bool run(const Term *t1, const Term *t2)
		{
			if (t1->key.elements_equal_to(t2->key))
            {
				PIRANHA_ASSERT(t1->cf.length() == 1 && t2->cf.length() == 1);

				return key_revlex_comparison_impl<typename Term::CfType::TermType>::run(&(*t1->cf.begin()), &(*t2->cf.begin()));

			} else
            {
				return t1->key.revlex_comparison(t2->key);
			}
		}
	};


	// Default value handler class. Suitable for use with POD integral types.
	template <class T>
	struct cm_value_handler
	{
		// Make really sure we use this only with integral types.
		PIRANHA_STATIC_CHECK(boost::is_integral<T>::value,"");
		// Assign to the minmax vector the values in the array key.
		template <class Key>
		void assign(std::vector<boost::numeric::interval<T> > &minmax, const Key &key)
		{
			PIRANHA_STATIC_CHECK((boost::is_same<T, typename Key::value_type>::value), "");

			typedef typename Key::size_type size_type;
			const size_type size = key.size();
			
            PIRANHA_ASSERT(size <= minmax.size());

			for (size_type i = 0; i < size; ++i) {
				minmax[i] = key[i];
			}
		}


		// Test the values in array key and, if they sit outside the corresponding minmax intervals,
		// update the intervals to include them.
		template <class Key>
		void test(std::vector<boost::numeric::interval<T> > &minmax, const Key &key)
		{
			PIRANHA_STATIC_CHECK((boost::is_same<T,typename Key::value_type>::value),"");
			typedef typename Key::size_type size_type;
			const size_type size = key.size();
			PIRANHA_ASSERT(size <= minmax.size());
			for (size_type i = 0; i < size; ++i) {
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
		// NOTE: this is used in coding, so it is probably harmless performance-wise to leave here the safe numeric_cast.
		max_fast_int to_max_fast_int(const T &x) const
		{
			return boost::numeric_cast<max_fast_int>(x);
		}


		// Post-decode hook.
		// NOTE: here we use implicit conversion because we know from the global representation analysis
		// that we won't have out-of-range conversions.
		void post_decode(T &value, const max_fast_int &decoded) const
		{
			value = decoded;
		}


		// Comparison operator is always true by default.
		bool operator==(const cm_value_handler &) const
		{
			return true;
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
			PIRANHA_STATIC_CHECK((boost::is_same<mp_rational,typename Key::value_type>::value),"");
			typedef typename Key::size_type size_type;
			const size_type size = key.size();
			
            PIRANHA_ASSERT(size <= minmax.size());
			
            for (size_type i = 0; i < size; ++i)
            {
				compute_new_lcd_and_update(minmax,key[i].get_den());
				// Assign current value.
				minmax[i] = key[i];
				// Assign tmp value.
				m_tmp = m_lcd;
				// Scale interval down to the lcd.
				minmax[i] *= m_tmp;
				PIRANHA_ASSERT(minmax[i].lower().get_den() == 1 && minmax[i].upper().get_den() == 1);
			}
		}
		template <class Key>
		void test(std::vector<boost::numeric::interval<mp_rational> > &minmax, const Key &key)
		{
			PIRANHA_STATIC_CHECK((boost::is_same<mp_rational,typename Key::value_type>::value),"");
			typedef typename Key::size_type size_type;
			const size_type size = key.size();
			PIRANHA_ASSERT(size <= minmax.size());
			for (size_type i = 0; i < size; ++i) {
				// Compute the new lcd and update the current minmax_vector.
				compute_new_lcd_and_update(minmax,key[i].get_den());
				// Assign tmp value.
				m_tmp = key[i];
				// Scale interval down to the new lcd.
				m_tmp *= m_lcd;
				PIRANHA_ASSERT(m_tmp.get_den() == 1);
				if (m_tmp < minmax[i].lower()) {
					minmax[i].assign(m_tmp,minmax[i].upper());
				} else if (m_tmp > minmax[i].upper()) {
					minmax[i].assign(minmax[i].lower(),m_tmp);
				}
				PIRANHA_ASSERT(minmax[i].lower().get_den() == 1 && minmax[i].upper().get_den() == 1);
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
		// Post-decode hook.
		void post_decode(mp_rational &value, const max_fast_int &decoded) const
		{
			value = boost::lexical_cast<mp_rational>(decoded);
			value /= m_lcd;
		}
		// Comparison will test the lcd.
		bool operator==(const cm_value_handler &other) const
		{
			return (m_lcd == other.m_lcd);
		}
		// Store current lcd into m_old_lcd, compute lcd between current lcd and input argument value,
		// store it into m_lcd and update the minmax vector to take into account the new lcd if needed.
		void compute_new_lcd_and_update(std::vector<boost::numeric::interval<mp_rational> > &minmax, const mp_integer &value)
		{
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
					PIRANHA_ASSERT(minmax[i].lower().get_den() == 1 && minmax[i].upper().get_den() == 1);
				}
			}
		}
		// Lowest common denominator.
		mp_integer	m_lcd;
		// Old value of the lcd.
		mp_integer	m_old_lcd;
		// Useful temporary value.
		mp_rational	m_tmp;
	};

	// Define a type to hold the min/max values of array keys in series.
	template <class Series, int N>
	struct cm_tuple_impl {
		PIRANHA_STATIC_CHECK(N > 0,"");
		// minmax type, to be used for limits of input series.
		typedef typename Series::TermType::key_type::value_type value_type;
		typedef boost::tuples::cons<std::vector<boost::numeric::interval<value_type> >,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_minmax> type_minmax;
		// mp_integer minmax type, to be used when calculating global limits.
		typedef boost::tuples::cons<std::vector<boost::numeric::interval<mp_integer> >,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_mp_minmax> type_mp_minmax;
		// max_fast_int minmax type.
		typedef boost::tuples::cons<std::vector<boost::numeric::interval<max_fast_int> >,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_max_fast_int_minmax> type_max_fast_int_minmax;
		// value_handler tuple.
		typedef boost::tuples::cons<cm_value_handler<value_type>,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_value_handler> type_value_handler;
		// Multi-precision coding tuple.
		typedef boost::tuples::cons<std::vector<mp_integer>,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_mp_coding_tuple> type_mp_coding_tuple;
		// Coding tuple.
		typedef boost::tuples::cons<std::vector<max_fast_int>,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_coding_tuple> type_coding_tuple;
		// Decoding tuple: vectors of pairs of max_fast_ints, num and den resp. of the decoding formula.
		typedef boost::tuples::cons<std::vector<std::pair<max_fast_int,max_fast_int> >,
			typename cm_tuple_impl<typename Series::TermType::CfType, N - 1>::type_decoding_tuple> type_decoding_tuple;
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
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_minmax type_minmax;
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_mp_minmax type_mp_minmax;
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_max_fast_int_minmax type_max_fast_int_minmax;
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_value_handler type_value_handler;
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_mp_coding_tuple type_mp_coding_tuple;
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_coding_tuple type_coding_tuple;
		typedef typename cm_tuple_impl<Series, Series::echelonLevel + 1>::type_decoding_tuple type_decoding_tuple;
	};

	template <class Series, class Tuple>
	struct cm_init_vector_tuples_impl {
		template <class ArgsTuple>
		static void run(const ArgsTuple &argsTuple, Tuple &t)
		{
			typedef typename Tuple::head_type::size_type size_type;
			t.get_head().resize(boost::numeric_cast<size_type>(argsTuple.template get<Series::TermType::key_type::position>().size()));
			cm_init_vector_tuples_impl<typename Series::TermType::CfType,typename Tuple::tail_type>::run(argsTuple,t.get_tail());
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
	inline void cm_init_vector_tuple(Tuple &t, const ArgsTuple &argsTuple)
	{
		PIRANHA_STATIC_CHECK(boost::tuples::length<Tuple>::value == Series::echelonLevel + 1, "");

		cm_init_vector_tuples_impl<Series, Tuple>::run(argsTuple, t);
	}


	template <class MinMaxTuple>
	struct cm_minmax2 {
		template <class Series, class ValueHandlerTuple>
		static void run_init(const Series &s, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			PIRANHA_STATIC_CHECK(boost::tuples::length<MinMaxTuple>::value == boost::tuples::length<ValueHandlerTuple>::value, "");
			PIRANHA_ASSERT(s.length() == 1);
			// NOTE: here key size could be less than the size of the vector in this tuple position,
			// hence the importance of having the vector default-initialised to zero. This happens because
			// we allow multiplication by series with fewer arguments.
			vh_tuple.get_head().assign(minmax_tuple.get_head(), s.begin()->key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_init(s.begin()->cf, minmax_tuple.get_tail(), vh_tuple.get_tail());
		}

		template <class Series, class ValueHandlerTuple>
		static void run_test(const Series &s, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			PIRANHA_ASSERT(s.length() == 1);
			vh_tuple.get_head().test(minmax_tuple.get_head(), s.begin()->key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_test(s.begin()->cf,minmax_tuple.get_tail(),vh_tuple.get_tail());
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
			PIRANHA_STATIC_CHECK(boost::tuples::length<MinMaxTuple>::value == boost::tuples::length<ValueHandlerTuple>::value, "");
			vh_tuple.get_head().assign(minmax_tuple.get_head(), term.key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_init(term.cf, minmax_tuple.get_tail(), vh_tuple.get_tail());
		}

		// Test minmax tuple with term's array representation, and, if needed, update min and max to include it.
		template <class Term, class ValueHandlerTuple>
		static void run_test(const Term &term, MinMaxTuple &minmax_tuple, ValueHandlerTuple &vh_tuple)
		{
			vh_tuple.get_head().test(minmax_tuple.get_head(),term.key);
			cm_minmax2<typename MinMaxTuple::tail_type>::run_test(term.cf,minmax_tuple.get_tail(),vh_tuple.get_tail());
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
			PIRANHA_STATIC_CHECK(boost::tuples::length<OpTuple>::value == boost::tuples::length<MinMaxTuple>::value,"");
			PIRANHA_STATIC_CHECK(boost::tuples::length<OpTuple>::value == boost::tuples::length<MpMinMaxTuple>::value,"");
			PIRANHA_STATIC_CHECK(boost::tuples::length<OpTuple>::value == boost::tuples::length<ValueHandlerTuple>::value,"");
			typedef typename MinMaxTuple::head_type::size_type size_type;
			PIRANHA_ASSERT(minmax1.get_head().size() == minmax2.get_head().size() && 
				global_minmax.get_head().size() == minmax2.get_head().size());
			// Harmonize the value handler tuples.
			ValueHandlerTuple::head_type::harmonize(vh1.get_head(),minmax1.get_head(),vh2.get_head(),minmax2.get_head());
			boost::numeric::interval<mp_integer> tmp1, tmp2;
			const size_type size = minmax1.get_head().size();
			for (size_type i = 0; i < size; ++i) {
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
						boost::numeric::hull(
							global_minmax.get_head()[i],
							tmp1 - tmp2
						),-tmp2
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
		static void run(mp_integer *prev_value, MpCt &mp_ct, const MpMinMaxTuple &mp_gt)
		{
			typedef typename MpCt::head_type::size_type size_type;
			// Assume same sizes.
			const size_type size = mp_ct.get_head().size();
			PIRANHA_ASSERT(size == mp_gt.get_head().size());
			for (size_type i = 0; i < size; ++i) {
				// Assign to the current value the previous one.
				mp_ct.get_head()[i] = *prev_value;
				// Multiply previous value by current minmax interval's width + 1.
				mp_integer width(boost::numeric::width(mp_gt.get_head()[i]));
				width += 1;
				*prev_value *= width;
			}
			cm_mp_ct_impl<typename MpCt::tail_type>::run(prev_value,mp_ct.get_tail(),mp_gt.get_tail());
		}
	};

	template <>
	struct cm_mp_ct_impl<boost::tuples::null_type> {
		static void run(mp_integer *, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

	// Compute coding tuple from tuple of multi-precision minmax vectors.
	template <class MpCt, class MpMinMaxTuple>
	void compute_mp_coding_tuple(MpCt &mp_ct, const MpMinMaxTuple &mp_gt)
	{
		PIRANHA_STATIC_CHECK(boost::tuples::length<MpCt>::value == boost::tuples::length<MpMinMaxTuple>::value,"");
		mp_integer prev_value(1);
		cm_mp_ct_impl<MpCt>::run(&prev_value,mp_ct,mp_gt);
	}

	template <class Tuple1, class Tuple2>
	struct tuple_vector_dot_impl {
		PIRANHA_STATIC_CHECK(boost::tuples::length<Tuple1>::value == boost::tuples::length<Tuple2>::value,"");
		typedef typename Tuple1::head_type::value_type value_type;
		typedef typename Tuple1::head_type::size_type size_type;
		static void run(const Tuple1 &t1, const Tuple2 &t2, value_type &retval)
		{
			const size_type size = t1.get_head().size();
			PIRANHA_ASSERT(size == t2.get_head().size());
			for (size_type i = 0; i < size; ++i) {
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
			const size_type size = mp_vector.size();
			PIRANHA_ASSERT(size == fast_vector.size());
			for (size_type i = 0; i < size; ++i) {
				fast_vector[i] = boost::lexical_cast<max_fast_int>(mp_vector[i]);
			}
		}
		static void run_impl(const std::vector<boost::numeric::interval<mp_integer> > &mp_vector,
			std::vector<boost::numeric::interval<max_fast_int> > &fast_vector)
		{
			typedef std::vector<boost::numeric::interval<mp_integer> >::size_type size_type;
			const size_type size = mp_vector.size();
			PIRANHA_ASSERT(size == fast_vector.size());
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

	template <class OpTuple, bool SubRequested>
	struct cm_code_impl2 {
		template <class CodingTuple, class Cf, class VhTuple>
		static void run(const CodingTuple &ct, const Cf &cf, const VhTuple &vh_tuple, max_fast_int &retval1, max_fast_int &retval2)
		{
			PIRANHA_STATIC_CHECK((boost::is_same<typename CodingTuple::head_type::value_type,max_fast_int>::value),"");
			PIRANHA_ASSERT(cf.length() == 1);

			typedef typename Cf::TermType::key_type::size_type size_type;
			const size_type size = cf.begin()->key.size();
			PIRANHA_ASSERT(size <= ct.get_head().size());
			max_fast_int tmp = 0;
			
            for (size_type i = 0; i < size; ++i)
            {
				tmp = ct.get_head()[i] * vh_tuple.get_head().to_max_fast_int(cf.begin()->key[i]);
				retval1 += tmp;

				if (!OpTuple::head_type::value)
                {
					retval2 -= tmp;
				} else if (SubRequested)
                {
					retval2 += tmp;
				}
			}
			cm_code_impl2<typename OpTuple::tail_type,SubRequested>::run(ct.get_tail(),cf.begin()->cf,vh_tuple.get_tail(),retval1,retval2);
		}
	};


	template <bool SubRequested>
	struct cm_code_impl2<boost::tuples::null_type,SubRequested> {
		template <class Cf>
		static void run(const boost::tuples::null_type &, const Cf &, const boost::tuples::null_type &, const max_fast_int &, const max_fast_int &) {}
	};


	template <class OpTuple>
	struct cm_code_impl1 {
		template <class CodingTuple, class Term, class VhTuple>
		static void run(const CodingTuple &ct, const Term &term, const VhTuple &vh_tuple, max_fast_int &retval1, max_fast_int &retval2)
		{
			PIRANHA_STATIC_CHECK((boost::is_same<typename CodingTuple::head_type::value_type,max_fast_int>::value),"");
			PIRANHA_ASSERT(term.key.size() <= ct.get_head().size());

			static const bool sub_requested = op_has_sub<OpTuple>::value;
			typedef typename Term::key_type::size_type size_type;
			max_fast_int tmp = 0;
			// NOTE: again the assumption that the sizes of vector and key are compatible. Need to sort this out...
			for (size_type i = 0; i < term.key.size(); ++i)
            {
				tmp = ct.get_head()[i] * vh_tuple.get_head().to_max_fast_int(term.key[i]);
				retval1 += tmp;
				if (!OpTuple::head_type::value) {
					retval2 -= tmp;
				} else if (sub_requested) {
					retval2 += tmp;
				}
			}
			cm_code_impl2<typename OpTuple::tail_type,sub_requested>::run(ct.get_tail(),term.cf,vh_tuple.get_tail(),retval1,retval2);
		}
	};

	// Code term using provided operations tuple type, coding tuple, value handler tuple and appending result to vector of codes v_codes.
	template <class OpTuple, class CodingTuple, class Term, class VhTuple>
	inline void cm_code(const CodingTuple &ct, const Term &term, const VhTuple &vh_tuple, max_fast_int &retval1, max_fast_int &retval2)
	{
		PIRANHA_STATIC_CHECK(boost::tuples::length<OpTuple>::value == boost::tuples::length<CodingTuple>::value,"");
		PIRANHA_STATIC_CHECK(boost::tuples::length<OpTuple>::value == boost::tuples::length<VhTuple>::value,"");
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
		PIRANHA_STATIC_CHECK(boost::tuples::length<DecodingTuple>::value == boost::tuples::length<MinMaxTuple>::value,"");
		max_fast_int prev_range = 1;
		cm_build_decoding_tuple_impl<DecodingTuple>::run(&prev_range,dt,minmax);
	}

	// This struct deals with coefficient series.
	template <class DecodingTuple>
	struct cm_decode_impl2 {
		template <class FinalCf, class MinMaxTuple, class Cf, class VhTuple, class ArgsTuple>
		static void run(const FinalCf &final_cf, const DecodingTuple &dt, const MinMaxTuple &gr, Cf &cf,
			const VhTuple &vh_tuple, const max_fast_int &code, const ArgsTuple &argsTuple)
		{
			PIRANHA_STATIC_CHECK((boost::is_base_of<BaseSeriesTag,Cf>::value),"");
			typedef typename Cf::TermType term_type;
			typedef typename term_type::key_type::size_type size_type;
			
            PIRANHA_ASSERT(dt.get_head().size() == argsTuple.template get<term_type::key_type::position>().size());
			PIRANHA_ASSERT(dt.get_head().size() == gr.get_head().size());

			const size_type size = dt.get_head().size();
			// This term is going to be inserted into the coefficient series.
			term_type term;
			term.key.resize(size);
			
            for (size_type i = 0; i < size; ++i)
            {
				const max_fast_int tmp = (code % dt.get_head()[i].first) / dt.get_head()[i].second +
					gr.get_head()[i].lower();
				vh_tuple.get_head().post_decode(term.key[i],tmp);
			}

			// Next recursion will operate on the term-to-be-inserted's coefficient.
			cm_decode_impl2<typename DecodingTuple::tail_type>::run(final_cf,dt.get_tail(),gr.get_tail(),term.cf,vh_tuple.get_tail(),code,argsTuple);
			// Before insering, let's make sure to clear up the contents of the coefficient series.
			cf.clearTerms();
			// Insert the result.
			cf.insert(term,argsTuple);
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
			const max_fast_int &code, const ArgsTuple &argsTuple)
		{
			PIRANHA_ASSERT(code >= 0);
			typedef typename Term::key_type::size_type size_type;
			// Make sure the sizes of the current tuples' elements are consistent.
			// (De)coding tuple sizes should always be the same as args tuple's.
			PIRANHA_ASSERT(dt.get_head().size() == argsTuple.template get<Term::key_type::position>().size());
			PIRANHA_ASSERT(dt.get_head().size() == gr.get_head().size());
			const size_type size = dt.get_head().size();
			term.key.resize(size);
			for (size_type i = 0; i < size; ++i) {
				const max_fast_int tmp = (code % dt.get_head()[i].first) / dt.get_head()[i].second +
					gr.get_head()[i].lower();
				vh_tuple.get_head().post_decode(term.key[i],tmp);
			}
			cm_decode_impl2<typename DecodingTuple::tail_type>::run(final_cf,dt.get_tail(),gr.get_tail(),term.cf,
				vh_tuple.get_tail(),code,argsTuple);
		}
	};

	// Decode given code into term. final_cf is the last coefficient in the echelon hierarchy, from outwards to inwards.
	template <class FinalCf, class DecodingTuple, class MinMaxTuple, class Term, class VhTuple, class ArgsTuple>
	inline void cm_decode(const FinalCf &final_cf, const DecodingTuple &dt, const MinMaxTuple &gr, Term &term, const VhTuple &vh_tuple,
		const max_fast_int &code, const max_fast_int &min_code, const ArgsTuple &argsTuple)
	{
		PIRANHA_STATIC_CHECK(boost::tuples::length<DecodingTuple>::value == boost::tuples::length<VhTuple>::value,"");
		PIRANHA_STATIC_CHECK(boost::tuples::length<DecodingTuple>::value == boost::tuples::length<MinMaxTuple>::value,"");
		cm_decode_impl1<DecodingTuple>::run(final_cf,dt,gr,term,vh_tuple,code - min_code,argsTuple);
	}
}

#endif
