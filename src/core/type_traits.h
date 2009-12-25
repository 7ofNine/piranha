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
#include <complex>

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

	/// is_lightweight type trait specialization for complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class.
	 */
	template <class T>
	struct is_lightweight<std::complex<T> >: is_lightweight<T>::type {};

	/// Default type trait for exact computations.
	/**
	 * Used to determine whether computations with the class are exact (e.g., integer and
	 * rational types) or not (e.g., floating-point). Defaults to false.
	 */
	template <class>
	struct is_exact: boost::false_type {};

	/// is_exact type trait specialization for complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class.
	 */
	template <class T>
	struct is_exact<std::complex<T> >: is_exact<T>::type {};

	/// Default type trait for trigonometric classes.
	/**
	 * Used to determine whether the class can represent basic trigonometric functions
	 * (sin and cos) exactly. Defaults to false.
	 */
	template <class>
	struct is_trig_exact: boost::false_type {};

	/// is_trig_exact type trait specialization for complex types.
	/**
	 * Inherits the type trait from the value_type of the complex class.
	 */
	template <class T>
	struct is_trig_exact<std::complex<T> >: is_trig_exact<T>::type {};

	// Recursively examine the echelon hierarchy of a series: if TypeTrait is true for all elements
	// of the hierarchy, then value is also true, otherwise it is false.
	template<class Series, int N, template <class> class TypeTrait, bool IsAnd>
	struct series_trait_impl {
		static const bool value =
			IsAnd ?
				(series_trait_impl<typename Series::term_type::cf_type,N - 1,TypeTrait,IsAnd>::value &&
				TypeTrait<typename Series::term_type::key_type>::value)
			:
				(series_trait_impl<typename Series::term_type::cf_type,N - 1,TypeTrait,IsAnd>::value ||
				TypeTrait<typename Series::term_type::key_type>::value)
			;
	};

	template<class Cf, template <class> class TypeTrait, bool IsAnd>
	struct series_trait_impl<Cf,0,TypeTrait,IsAnd> {
		static const bool value = TypeTrait<Cf>::value;
	};

	/// Type-trait helper for series.
	/**
	 * This struct defines a static const boolean value in the following way: if all (or any) of the elements of the series hierarchy
	 * have satisfied the type trait TypeTrait, then value is true, otherwise value is false.
	 * The choice between logical OR and logical AND is dictated by the boolean flag IsAnd.
	 */
	template <class Series, template <class> class TypeTrait, bool IsAnd>
	struct series_trait {
		/// Static boolean value.
		static const bool value = series_trait_impl<Series,Series::echelon_level + 1,TypeTrait,IsAnd>::value;
	};
}

#endif
