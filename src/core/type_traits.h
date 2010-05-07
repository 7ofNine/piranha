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

#ifndef PIRANHA_TYPE_TRAITS_H
#define PIRANHA_TYPE_TRAITS_H

#include <boost/type_traits/integral_constant.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>

#include "base_classes/base_series_tag.h"
#include "config.h"

namespace piranha
{
	// Fwd declaration.
	template <class Cf, class Key> class tetd_helper;

	template <class Term>
	class term_eval_type_determiner
	{
		public:
			typedef typename tetd_helper<typename Term::cf_type::eval_type, typename Term::key_type::eval_type>::type type;
	};

	template <class CfEval, class KeyEval>
	class tetd_helper
	{
		public:
			typedef double type;
	};

	template <>
	class tetd_helper<double,std::complex<double> >
	{
		public:
			typedef std::complex<double> type;
	};

	template <>
	class tetd_helper<std::complex<double>,double>
	{
		public:
			typedef std::complex<double> type;
	};

	template <>
	class tetd_helper<std::complex<double>,std::complex<double> >
	{
		public:
			typedef std::complex<double> type;
	};

	/// Default lightweight type trait.
	/**
	 * Used to determine if it is worth to cache objects in contiguous memory areas during
	 * performance-critical sections. Defaults to false.
	 */
	template <class>
	struct is_lightweight: boost::false_type {};

	/// is_lightweight type trait specialisation for complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class.
	 */
	template <class T>
	struct is_lightweight<std::complex<T> >: is_lightweight<T>::type {};

	/// Default type trait for exact ring operations.
	/**
	 * Used to determine whether ring operations (add, sub and mult) with the class are exact (e.g., integer and
	 * rational types) or not (e.g., floating-point types). Defaults to false.
	 */
	template <class, class Enable = void>
	struct is_ring_exact: boost::false_type {};

	/// is_ring_exact type trait specialisation for non-series complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class. Specialisation is disabled
	 * if T is not a series type (in that case, the series type trait specialisation will be used).
	 */
	template <class T>
	struct is_ring_exact<T,typename boost::enable_if_c<boost::is_complex<T>::value && !boost::is_base_of<base_series_tag,T>::value>::type>:
		is_ring_exact<typename T::value_type>
	{};

	/// is_ring_exact type trait specialisation for series.
	/**
	 * Will be true if coefficient and key are ring exact.
	 */
	template <class T>
	struct is_ring_exact<T,typename boost::enable_if<boost::is_base_of<base_series_tag,T> >::type>
	{
		static const bool value = is_ring_exact<typename T::term_type::cf_type>::value && is_ring_exact<typename T::term_type::key_type>::value;
	};

	/// Default type trait for trigonometric classes.
	/**
	 * Used to determine whether the class can represent basic trigonometric functions
	 * (sin and cos) exactly. Defaults to false.
	 */
	template <class, class Enable = void>
	struct is_trig_exact: boost::false_type {};

	/// is_trig_exact type trait specialisation for non-series complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class. Specialisation is disabled
	 * if T is not a series type (in that case, the series type trait specialisation will be used).
	 */
	template <class T>
	struct is_trig_exact<T,typename boost::enable_if_c<boost::is_complex<T>::value && !boost::is_base_of<base_series_tag,T>::value>::type>:
		is_trig_exact<typename T::value_type>
	{};

	/// is_trig_exact type trait specialisation for series.
	/**
	 * Will be true if either coefficient or key are trig exact.
	 */
	template <class T>
	struct is_trig_exact<T,typename boost::enable_if<boost::is_base_of<base_series_tag,T> >::type>
	{
		static const bool value = is_trig_exact<typename T::term_type::cf_type>::value || is_trig_exact<typename T::term_type::key_type>::value;
	};

	/// Default type trait for classes dividable by int.
	/**
	 * Used to determine whether the class can be divided exactly by an int
	 * (e.g., rationals). Defaults to false.
	 */
	template <class, class Enable = void>
	struct is_divint_exact: boost::false_type {};

	/// is_divint_exact type trait specialisation for non-series complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class. Specialisation is disabled
	 * if T is not a series type (in that case, the series type trait specialisation will be used).
	 */
	template <class T>
	struct is_divint_exact<T,typename boost::enable_if_c<boost::is_complex<T>::value && !boost::is_base_of<base_series_tag,T>::value>::type>:
		is_divint_exact<typename T::value_type>
	{};

	/// is_divint_exact type trait specialisation for series.
	/**
	 * Will be true if either coefficient or key are divint exact.
	 */
	template <class T>
	struct is_divint_exact<T,typename boost::enable_if<boost::is_base_of<base_series_tag,T> >::type>
	{
		static const bool value = is_divint_exact<typename T::term_type::cf_type>::value || is_divint_exact<typename T::term_type::key_type>::value;
	};

	template <class CfSeries, class Enable = void>
	struct final_cf_impl
	{
		typedef typename final_cf_impl<typename CfSeries::term_type::cf_type>::type type;
	};

	template <class Cf>
	struct final_cf_impl<Cf,typename boost::enable_if_c<!boost::is_base_of<base_series_tag,Cf>::value>::type>
	{
		typedef Cf type;
	};

	/// Final coefficient type.
	/**
	 * Coefficient type at the end of the echelon recursion.
	 */
	template <class Series>
	struct final_cf
	{
		p_static_check((boost::is_base_of<base_series_tag,Series>::value),"Cannot determine final coefficient of a non-series type.");
		typedef typename final_cf_impl<typename Series::term_type::cf_type>::type type;
	};
}

#endif
