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


#include "base_classes/base_series_tag.h"
#include "config.h"
#include "mp.h"
#include "base_classes/numerical_container_tag.h"

#include <boost/type_traits/integral_constant.hpp>
#include <boost/type_traits/is_complex.hpp>

#include <complex>
#include <type_traits>

namespace piranha
{
	// Fwd declaration.
	template <class Cf, class Key> class TermEvalTypeDeterminerHelper;

    // templates for determining type  to which a term evaluates when a time value is used
	template <class Term>
	class TermEvalTypeDeterminer
	{
		public:
			typedef typename TermEvalTypeDeterminerHelper<typename Term::CfType::EvalType, typename Term::KeyType::EvalType>::Type Type;
	};

	template <class CfEval, class KeyEval>
	class TermEvalTypeDeterminerHelper
	{
		public:
			typedef double Type;
	};

	template <>
	class TermEvalTypeDeterminerHelper<double, std::complex<double> >
	{
		public:
			typedef std::complex<double> Type;
	};

	template <>
	class TermEvalTypeDeterminerHelper<std::complex<double> ,double>
	{
		public:
			typedef std::complex<double> Type;
	};

	template <>
	class TermEvalTypeDeterminerHelper<std::complex<double>, std::complex<double> >
	{
		public:
			typedef std::complex<double> Type;
	};


	// Piranha Series type
	template <typename T>
	concept  PiranhaSeries = std::is_base_of_v<BaseSeriesTag, T>;

	//template <typename T, typename Derived>
	//concept PiranhaSimilarSeries = requires {
	//	std::conjunction < std::is_base_of<BaseSeriesTag, T>, std::disjunction <std::is_same<Derived, T>, std::is_same<Derived, std::complex<T> > >};
	////std::is_base_of_v<BaseSeriesTag, T> && (std::is_same_v<Derived, T> || std::is_same_v<Derived, std::complex<T> >)


	/// Default type trait for exact ring operations.
	/**
	 * Used to determine whether ring operations (add, sub and mult) with the class are exact (e.g., integer and
	 * rational types) or not (e.g., floating-point types). Defaults to false.
	 */
	template <typename>
	struct is_ring_exact: boost::false_type {};

	/// is_ring_exact type trait specialisation for non-series complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class. Specialisation is disabled
	 * if T is not a series type (in that case, the series type trait specialisation will be used).
	 */
	template <typename T> requires boost::is_complex<T>::value && (!PiranhaSeries<T>)
	struct is_ring_exact<T>:
		is_ring_exact<typename T::value_type>
	{};

	/// is_ring_exact type trait specialisation for series.
	/**
	 * Will be true if coefficient and key are ring exact.
	 */
	template <PiranhaSeries T>
	struct is_ring_exact<T> 
	{
		static const bool value = is_ring_exact<typename T::TermType::CfType>::value && is_ring_exact<typename T::TermType::KeyType>::value;
	};


	/// Default type trait for trigonometric classes.
	/**
	 * Used to determine whether the class can represent basic trigonometric functions
	 * (sin and cos) exactly. Defaults to false.
	 */
	template <typename>
	struct is_trig_exact: boost::false_type {};

	/// is_trig_exact type trait specialisation for non-series complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class. Specialisation is disabled
	 * if T is not a series type (in that case, the series type trait specialisation will be used).
	 */
	template <typename T> requires boost::is_complex<T>::value && (!PiranhaSeries<T>)
	struct is_trig_exact<T>:
		is_trig_exact<typename T::value_type>
	{};

	/// is_trig_exact type trait specialisation for series.
	/**
	 * Will be true if either coefficient or key are trig exact.
	 */
	template <PiranhaSeries T>
	struct is_trig_exact<T>
	{
		static const bool value = is_trig_exact<typename T::TermType::CfType>::value || is_trig_exact<typename T::TermType::KeyType>::value;
	};

	/// Default type trait for classes dividable by int.
	/**
	 * Used to determine whether the class can be divided exactly by an int
	 * (e.g., rationals). Defaults to false.
	 */
	template <typename>
	struct is_divint_exact: boost::false_type {};

	/// is_divint_exact type trait specialisation for non-series complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class. Specialisation is disabled
	 * if T is not a series type (in that case, the series type trait specialisation will be used).
	 */
	template <typename T> requires boost::is_complex<T>::value && (!PiranhaSeries<T>)
	struct is_divint_exact<T>:
		is_divint_exact<typename T::value_type>
	{};

	/// is_divint_exact type trait specialisation for series.
	/**
	 * Will be true if either coefficient or key are divint exact.
	 */
	template <PiranhaSeries T>
	struct is_divint_exact<T>
	{
		static const bool value = is_divint_exact<typename T::TermType::CfType>::value || is_divint_exact<typename T::TermType::KeyType>::value;
	};

	/// Default type trait for classes which can represent rational exponents.
	/**
	 * Used to determine whether the class can represent symbolic quantities raised to rational exponents.
	 * Defaults to false.
	 */
	template <typename>
	struct is_rational_exponent: boost::false_type {};

	/// is_rational_exponent type trait specialisation for series.
	/**
	 * Will be true if series has a degree_type typedef which is rational.
	 */
	template <PiranhaSeries T> requires std::is_same_v<typename T::degree_type, mp_rational>
	struct is_rational_exponent<T>:
	boost::true_type {};


    //
    // traits for final coefficient type of an echeloned series
	//
    template <typename T> 
	struct FinalCfImplementation
	{
		using Type = typename FinalCfImplementation<typename T::TermType::CfType>::Type;
	};


	template <typename T> requires (!PiranhaSeries<T>)
	struct FinalCfImplementation<T> 
	{
		using Type = T;
	};




	/// Final coefficient type.
	/**
	 * Coefficient type at the end of the echelon recursion.
	 */
	template <PiranhaSeries T>
	struct FinalCf
	{
		using Type = typename FinalCfImplementation<typename T::TermType::CfType>::Type;
	};

	template <typename T>
	concept PiranhaNumerical = std::is_base_of_v<numerical_container_tag, T>;
}

#endif
