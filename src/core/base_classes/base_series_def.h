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

#ifndef PIRANHA_BASE_SERIES_DEF_H
#define PIRANHA_BASE_SERIES_DEF_H

#include <boost/functional/hash.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/unordered_set.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_pointer.hpp>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../memory.h"
#include "../mp.h"
#include "../Psym.h"
#include "../type_traits.h"
#include "base_series_mp.h"
#include "base_series_tag.h"

// Template parameters list for piranha::BaseSeries (declaration form).
#define __PIRANHA_BASE_SERIES_TP_DECL class Term, char Separator, class Allocator, class Derived
// Template parameters list for piranha::BaseSeries (implementation form).
#define __PIRANHA_BASE_SERIES_TP Term, Separator, Allocator, Derived

namespace piranha
{
	// Implementation of echelon level determination.
	template <class CfSeries, class Enable = void>
	struct EchelonLevelImpl
	{
		static const int value = EchelonLevelImpl<typename CfSeries::TermType::CfType>::value + 1;
	};


	template <class FinalCf>
	struct EchelonLevelImpl<FinalCf, typename boost::enable_if_c<!boost::is_base_of<BaseSeriesTag, FinalCf>::value>::type>
	{
		static const int value = 0;
	};


	// Struct to define the container used in series - like a template typedef.
	template <class Term>
	struct SeriesContainer
	{
		typedef boost::unordered_set<Term, boost::hash<Term>, std::equal_to<Term>, CountingAllocator<Term, std::allocator<Term> > > Type;
	};


	// These accessors are used in generic code that must work on both plain series (i.e., iterators) and sorted representations
	// of series as returned by get_sorted_series (i.e., pointers to pointers of terms).
	template <class Iterator, class Enable = void>
	struct FromIterator
	{
		static const typename Iterator::value_type *get(const Iterator &it)
		{
			return &(*it);
		}
	};


	template <class Iterator>
	struct FromIterator<Iterator, typename boost::enable_if<boost::is_pointer<typename Iterator::value_type> >::type>
	{
		static typename Iterator::value_type get(const Iterator &it)
		{
			return *it;
		}
	};



	/// Base series class.
	/**
	 * This class provides the basic representation of a series as a collection of terms stored into a hash set. The class is intended
	 * to be inherited together with (at least) either piranha::NamedSeries (for a top-level series) or piranha::cf_series (for a coefficient series).
	 *
	 * The methods in this class provide the lowest level of series manipulation, and allow to operate directly on the individual terms of the
	 * series.
	 *
	 * @author Francesco Biscani (bluescarni@gmail.com)
	 */

    //
    // for series multiplication it requires to be used in parallel with series_multiplcation which in turn requires
    // a multiplier and a truncator
    //
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	class BaseSeries: BaseSeriesTag
	{
			// Befriend meta-programming classes.
			template <class, class>
			friend struct SeriesFromCfImpl;

			template <class, class>
			friend struct SeriesFromKeyImpl;

			template <int>
			friend class SeriesFlattener;

			template <class, class>
			friend struct BaseSeriesAddSelector;

			template <class, class>
			friend struct BaseSeriesSubtractSelector;

			template <class, class, class>
			friend struct BaseSeriesMultiplySelector;

			template <class, class>
			friend struct BaseSeriesEqualToSelector;

		public:

			/// Alias for term type.
			typedef Term TermType;

			/// Term container.
			/**
			 * The underlying term container is a plain boost::unordered set. Term types must specialise the boost::hash class, which
			 * will be used to provide the hash values for terms.
			 */
			typedef typename SeriesContainer<TermType>::Type ContainerType;

			/// Echelon level.
			static const int echelonLevel = EchelonLevelImpl<typename TermType::CfType>::value;

			/// Type resulting from series evaluation.
			/**
			 * Determined automatically at compile time. Will be double in case of series with scalar coefficients,
			 * std::complex<double> in case of series with complex coefficients.
			 */
			typedef typename TermEvalTypeDeterminer<Term>::Type EvalType;

			/// Term separator in textual series representation.
			static const char separator = Separator;

			/// Const iterator over the series' terms.
			typedef typename ContainerType::const_iterator const_iterator;

			/// Size type.
			typedef typename ContainerType::size_type size_type;

			/** @name Series properties. */
			//@{
			size_type length()       const;
			bool      empty()        const;
			bool      isSingleCf()   const;
			size_type atoms()        const;
			//@}


		//protected:
			/** @name Construction/destruction. */
			//@{
            //defined: base_series_io.h
			template <class Key, class ArgsTuple>
			static Derived baseSeriesFromKey(const Key &, const ArgsTuple &);

			template <class Cf, class ArgsTuple>
			static Derived baseSeriesFromCf(const Cf &, const ArgsTuple &);

			template <class Number, class ArgsTuple>
			static Derived baseSeriesFromNumber(const Number &, const ArgsTuple &);

			template <class ArgsTuple>
			void baseConstructFromPsym(const Psym &, const int, const ArgsTuple &);

			~BaseSeries();
			//@}

			/** @name Series manipulation. */
			//@{
			void eraseTerm(const const_iterator &);

			void clearTerms();

			template <bool, bool, class Term2, class ArgsTuple>
			void insert(const Term2 &, const ArgsTuple &);

			template <class Term2, class ArgsTuple>
			void insert(const Term2 &, const ArgsTuple &);

			template <class Iterator, class ArgsTuple>
			void insertRange(const Iterator &, const Iterator &, const ArgsTuple &);

			void baseSwap(Derived &);

            //
            // not sure what this is good for. requires series_multiplication and truncator
            // I don't think it is used anywhere. except to define split in cf_series and named_series
            // 
			template <class Series, class ArgsTuple>
			void baseSplit(std::vector<std::vector<Series> > &, const int n, const ArgsTuple &) const;

            //
            // don't understand how that one really works.
            // test with Echelon level > 0
            //
            // create a vector of the series term. Basically use distributive law.
            // (x+2y+z)*cos(a) -> vector(x*cos, 2y*cos(a), z*cos(a))
            //
			template <class ArgsTuple>
			std::vector<TermType> flattenTerms(const ArgsTuple &) const;

            //
            // layouts apply to named series and describe the relation of the arguments to the key
            //
            //
            // applies a new layout to the terms of the BaseSeries. LayoutTuple is tuple (number of echonLevel members)
            // of vectors of pairs (bool, int). 
            // In BaseSeries we iterate over the the terms of the BaseSeries term container. 
            // for the cf as well as the key of the term we iterate over the size of the vector in the tuple and create a new
            // container of this size. The size of the current BasSeries term cf/key container has to be smaller or equal
            // to the size of the echelonLevel corresponding vector of the LayoutTuple. If the pair in the vector is 
            // (true,int) int gives an index into the current container and this index will be transfered to
            // new container at the current index position. Consequentially the int must be samller (it is an index) than the size of the current 
            // BaseSeries container. 
            //
			template <class LayoutTuple, class ArgsTuple>
			void applyLayoutToTerms(LayoutTuple const &, ArgsTuple const &, Derived &) const;

            //
            // test if we can trim some arguments from the argsTuple and terms. They can be trimmed if 
            // return the result as trimFlags(tuple of vector<bool>, false: trim i.e. the symbol in the key of the term is not used
            // The trimFlags vector has to be allocated and sized before calling trimTestTerms and has to agree with the sizes of the
            // key(s) used for the terms.
            // Used in cf_series and named_series.
            //
			template <class TrimFlags>
			void trimTestTerms(TrimFlags &) const;

            //
            // trim the keys according to the trimFlags from the terms and return the updated series in retval
            // The change of the argumentsTuple was done earlier before the use of trimTerms and the return series is setup
            // accordingly.
            //
			template <class TrimFlags, class ArgsTuple>
			void trimTerms(TrimFlags const &, ArgsTuple const &, Derived &) const;

            //
            // To be clarified
            //
			template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries baseSub(const PosTuple &, SubCaches &, const ArgsTuple &) const;
            //@}

			/** @name Terms accessors. */
			//@{
            // 
            // find specified term in the BaseSeries term container
            //
            const_iterator findTerm(const TermType &) const; // definition: base_series_manip.h
			
            //
            // first element of container
            //
            const_iterator begin() const;

            //
            // end of container 
            //
			const_iterator end() const;
			//@}

			/** @name Base comparison methods. */
			//@{
			template <class T>
			bool baseEqualTo(const T &) const;

			//@}
			/** @name Base maths. */
			//@{
            //
            // calculate the norm of the series
            //
			template <class ArgsTuple>
			double baseNorm(const ArgsTuple &) const;

            //
            // evaluate series at value (a numeric value)
            // result (EvalType) can be double or complex<double> depending on series
            //
			template <class ArgsTuple>
			EvalType baseEval(const double, const ArgsTuple &) const;

			template <class T, class ArgsTuple>
			Derived &baseAdd(const T &, const ArgsTuple &);

			template <class T, class ArgsTuple>
			Derived &baseSubtract(const T &, const ArgsTuple &);

            //
            // requires series_multiplication to be usefull
            //
			template <class T, class ArgsTuple>
			Derived &baseMultBy(const T &, const ArgsTuple &);

			template <class T, class ArgsTuple>
			Derived &baseDivideBy(const T &, const ArgsTuple &);

			template <class ArgsTuple>
			Derived baseInvert(const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived basePow(const double &, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived basePow(const mp_rational &, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived naturalPower(const std::size_t, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived negativeIntegerPower(const int, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived realPower(const double, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived rationalPower(const mp_rational &, const ArgsTuple &) const;

			template <class ArgsTuple>
			Derived baseRoot(const int, const ArgsTuple &) const;

			template <class Series, class PosTuple, class ArgsTuple>
			static void basePartial(const Derived &, Series &, const PosTuple &, const ArgsTuple &);

			template <class PosTuple, class ArgsTuple>
			Derived basePartial(int, const PosTuple &, const ArgsTuple &) const;

			template <class PosTuple, class ArgsTuple>
			Derived basePartial(const PosTuple &, const ArgsTuple &) const;
			//@}

			/** @name Base output streaming methods. */
			//@{
			template <class ArgsTuple>
			void printTermsPlain(std::ostream &, const ArgsTuple &) const;

			template <class ArgsTuple>
			void printTermsTEX(std::ostream &, const ArgsTuple &) const;

			template <class ArgsTuple>
			void printTermsPretty(std::ostream &, const ArgsTuple &) const;
            //@}

			// Standard substitution functor. Will call sub() on coefficients and keys.
			class SubstitutionFunctor
            {
                public:

				template <class RetSeries, class Element, class PositionTuple, class SubsitutionCaches, class ArgsTuple>
				static RetSeries run(Element const &element, PositionTuple const &positionTuple, SubsitutionCaches &substitutionCaches, ArgsTuple const &argsTuple)
				{
					return element.template sub<RetSeries>(positionTuple, substitutionCaches, argsTuple);
				}
			};

		private:

			template <class T, class ArgsTuple>
			Derived &multiplyCoefficientsBy(const T &, const ArgsTuple &);

			template <class T, class ArgsTuple>
			Derived &divideCoefficientsBy(const T &, const ArgsTuple &);

			template <int, class T, class ArgsTuple>
			void multDivCoefficientsBy(const T &, const ArgsTuple &);

			template <bool, class Number, class ArgsTuple>
			Derived &mergeWithNumber(const Number &, const ArgsTuple &);

			template <bool, class Derived2, class ArgsTuple>
			Derived &mergeTerms(const Derived2 &, const ArgsTuple &);

			//
            // compare series on base level. T should be a BaseSeries or numerical
            //
            template <class T>
			bool genericSeriesComparison(const T &) const;

			template <class Number>
			bool genericNumericalComparison(const Number &) const;

			template <class Iterator, class ArgsTuple>
			void genericPrintTermsPretty(std::ostream &, const Iterator &, const Iterator &, const ArgsTuple &) const;

			template <class Iterator, class ArgsTuple>
			void genericPrintTermsTEX(std::ostream &, const Iterator &, const Iterator &, const ArgsTuple &) const;

			template <class Iterator, class Series, class ArgsTuple>
			void genericBaseSplit(std::vector<std::vector<Series> > &, const Iterator &, const Iterator &, const ArgsTuple &) const;

			template <bool, class ArgsTuple>
			void llInsert(const TermType &, const ArgsTuple &);

			template <bool, class ArgsTuple>
			void termInsertNew(const TermType &, const ArgsTuple &);

			template <class Number, class ArgsTuple>
			bool commonPowHandler(const Number &, Derived &retval, const ArgsTuple &) const;

		private:

			ContainerType container;
	};

#define E0_SERIES_TP_DECL class Cf, class Key, class Multiplier, class Truncator, class Allocator
#define E0_SERIES_TP Cf, Key ,Multiplier, Truncator, Allocator
#define E0_SERIES_TERM(TermName) TermName<Cf, Key, '|', Allocator>
#define E0_SERIES(SeriesName) SeriesName<E0_SERIES_TP>
#define E0_SERIES_BASE_ANCESTOR(TermName, SeriesName) piranha::BaseSeries<E0_SERIES_TERM(TermName), '\n', \
	Allocator, E0_SERIES(SeriesName) >

#define E1_SERIES_TP_DECL class Cf, class Key0, class Key1, \
						  class Mult0, class Mult1, class Trunc0, class Trunc1, class Allocator
#define E1_SERIES_TP      Cf, Key0, Key1, Mult0, Mult1, Trunc0, Trunc1, Allocator
#define E1_SERIES_COEFFICIENT(CfName)     CfName<Cf, Key0, Mult0, Trunc0, Allocator>
#define E1_SERIES(SeriesName)             SeriesName<E1_SERIES_TP>
#define E1_SERIES_TERM(TermName, CfName) TermName< CfName, Key1, '|', Allocator >
#define E1_SERIES_BASE_ANCESTOR(TermName, CfName, SeriesName) piranha::BaseSeries<TermName< \
	CfName, Key1, '|', Allocator>, \
	'\n', Allocator, SeriesName >
}

#endif
