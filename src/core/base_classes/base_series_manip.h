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
#include <iostream>

#include "../config.h"
#include "../exceptions.h"
#include "base_series_def.h"
#include "base_series_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Transform term.
	/**
	 * Non-specialized version, it will create a copy of the converted term.
     *
     * converts from type Term2 to type Term1 by constructing a new term with possible change for the arguments
     * it the types are the same  a copy is created (see spcialization below)
	 */
	template <class Term2, class Term1>
	class TermConverter
	{
		public:
			/// Constructor.
			template <class ArgsTuple>
			explicit TermConverter(const Term2 &c, const ArgsTuple &a): result(c, a) {}
			/// "Copy" of the converted term in terms of type Term1
			const Term1 result;
	};


	/// Specialized term converter.
	/**
	 * It will be invoked when the type to convert from is the same as the converted type. A reference
	 * to the convertee is stored inside the class.
	 */
	template <class Term2>
	class TermConverter<Term2, Term2>
	{
		public:
			/// Constructor.
			template <class ArgsTuple>
			explicit TermConverter(const Term2 &c, const ArgsTuple &): result(c) {}
			/// Reference to the converted term.
			const Term2 &result;
	};


	// High-level insertion function.
	//
	// This function is used to insert terms into a series. It requires that the number of arguments
	// of each element of the term is smaller or equal to the series',
	// otherwise an assertion fails and the program aborts. base_pseries::mergeArgs,
	// base_pseries::append_cf_args, base_pseries::append_trig_args, etc. can be used to add the needed arguments
	// to the series.
	//
	// This function performs some checks and then calls llInsert.
	//
	// CanonicalCheck: true: a canonicality check is done for the key of term. For an exponential key this means no effect.
	//                       For a trigonometric key this means that the first (?) argument in the key is not negative, which
	//						 might mean changing the signs of the coefficent and keys (sinus) of the term
	// Sign: true : add the term , for an insert just add the positive value of the term
	//	     false: subtract the term. For a new insert that means insert the negative value of the term.
	//       Warning: this is somewhat confusing!!
	//
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <bool CanonicalCheck, bool Sign, class Term2, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term2, const ArgsTuple &argsTuple)
	{
#ifdef DEBUG
        bool const TCanonical = CanonicalCheck;
        bool const TSign = Sign;
#endif
        //convert the inserted term to the type of the term (TermType) of the series it is being inserted to
		TermConverter<Term2, TermType> convertedTerm(term2, argsTuple);

		// Make sure the appropriate routines for the management of arguments have been called.
		PIRANHA_ASSERT(convertedTerm.result.cf.isInsertable(argsTuple) && convertedTerm.result.key.isInsertable(argsTuple));

		TermType *newTerm = nullptr;
		if (convertedTerm.result.cf.needsPadding(argsTuple) || convertedTerm.result.key.needsPadding(argsTuple)) 
		{
			newTerm = CountingInterface::allocate(allocator, 1, 0);
			CountingInterface::construct(allocator, newTerm, convertedTerm.result);
			newTerm->cf.padRight(argsTuple);
			newTerm->key.padRight(argsTuple);
		}

		if (CanonicalCheck) 
		{
			if (!convertedTerm.result.isCanonical(argsTuple)) 
			{
				if (newTerm == 0) 
				{
					newTerm = CountingInterface::allocate(allocator, 1, 0);
					CountingInterface::construct(allocator, newTerm, convertedTerm.result);
				}
				newTerm->canonicalise(argsTuple);
			}
		}

		const TermType *insertTerm = newTerm != nullptr ? newTerm : &convertedTerm.result;

		llInsert<Sign>(*insertTerm, argsTuple);
		
		if (newTerm != nullptr) 
		{
			CountingInterface::destroy(allocator, newTerm);
			CountingInterface::deallocate(allocator, newTerm, 1);
		}
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Term2, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term, const ArgsTuple &argsTuple)
	{
		insert<true, true>(term, argsTuple);
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::insertRange(const Iterator &begin, const Iterator &end, const ArgsTuple &argsTuple)
	{
		for (Iterator it = begin; it != end; ++it) 
		{
			insert(*it, argsTuple);
		}
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<PIRANHA_BASE_SERIES_TP>::findTerm(const TermType &term) const
	{
		return container.find(term);
	}


    //the low level insert routine. All the adjustments are out of the way. See the PIRANHA_ASSERT
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::llInsert(const TermType &term, const ArgsTuple &argsTuple)
	{
#ifdef DEBUG
        bool TSign = Sign;
#endif
		// TODO: think about moving this check higher in the stack of functions, we probably don't want to reach
		// _this_ point before checking for ignorability.
		if (term.cf.isIgnorable(argsTuple) || term.key.isIgnorable(argsTuple)) 
		{
			return;
		}
		PIRANHA_ASSERT(term.cf.isInsertable(argsTuple) && term.key.isInsertable(argsTuple) &&
			          !term.cf.needsPadding(argsTuple) && !term.key.needsPadding(argsTuple) && term.isCanonical(argsTuple));

		const_iterator it(findTerm(term));
		if (it == end()) 
		{
			// The term is NOT a duplicate, insert in the set.
			termInsertNew<Sign>(term, argsTuple);

		} else 
		{
			// The term is in the set, hence an existing term will be modified.
			// Add or subtract according to request.
			if (Sign) 
			{
				it->cf.add(term.cf, argsTuple);

			} else 
			{
				it->cf.subtract(term.cf, argsTuple);
			}

			// Check if the new resulting coefficient can be ignored. Beware there would have been an addition or subtraction.
			// in case of new insertion the term would already have been dumped at the beginning.
			if (it->cf.isIgnorable(argsTuple)) 
			{
				eraseTerm(it);
			}
		}
	}


	// Insert a new term into the series
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::termInsertNew(const TermType &term, const ArgsTuple &argsTuple)
	{
#ifdef DEBUG
        bool TSign = Sign;
#endif
        //insert the term into the container
		std::pair<const_iterator, bool> res(container.insert(term));

		PIRANHA_ASSERT(res.second); // check that insert worked

        //if not sign invert the coefficient
		//if Sign is true than it is like adding an element
		//if Sign is fale we have a subtraction. 
		// Inserting a newTerm should normaly handled like adding a new term.
		// sometimes we might handle it as inserting the negative value.
		if (!Sign) // was that ever the other way round?? what is actually intended ??
		{
			res.first->cf.invertSign(argsTuple);
		}
	}

    //remove term given by the iterator
	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::eraseTerm(const const_iterator &it)
	{
		container.erase(it);
	}


	/// Swap the terms with another series.
	/**
	 * All terms get swapped.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSwap(Derived &series2)
	{
		PIRANHA_ASSERT(derived_cast != &series2);
		container.swap(series2.container);  // swap the series terms by swapping the terms container
	}

    //
	// Apply an arguments layout to all terms and insert them into retval.
    // for layout see named_series_manip.h.
    // The control over the proper assigment of the positions to the arguments is outisde of this method somewhere in named series. Only there
    // layouts actually make sense.
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class LayoutTuple, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::applyLayoutToTerms(LayoutTuple const &layoutTuple, ArgsTuple const &argsTuple, Derived &retval) const
	{
		for (TermType term : *this)
		{
			term.cf.applyLayout(layoutTuple, argsTuple);
			term.key.applyLayout(layoutTuple, argsTuple);
			retval.insert(term, argsTuple);
		}
	}


    //
    // test if we can trim some arguments from the argsTuple and terms
    // return the result as trimFlags(tuple of vector<bool>, false: trim i.e. the symbol in the argumentTuple is not used
    //
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::trimTestTerms(TrimFlags &trimFlags) const
	{
		const const_iterator itEnd = end();
		for (const_iterator it = begin(); it != itEnd; ++it) 
		{
            // this is accumulative over the series
			it->cf.trimTest(trimFlags);
			it->key.trimTest(trimFlags);
		}
	}


    //
    // trim the arguments according to the trimFlags from the *this terms argumentTuple and return the modified series in retval
    // the trimmming of the argumentsTuple was done earlier
    //
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::trimTerms(TrimFlags const &trimFlags, ArgsTuple const &argsTuple, Derived &retval) const
	{
		const const_iterator itEnd = end();
		for (const_iterator it = begin(); it != itEnd; ++it) 
		{
			retval.insert(TermType(typename TermType::CfType(it->cf.trim(trimFlags, argsTuple)), typename TermType::KeyType(it->key.trim(trimFlags, argsTuple))), argsTuple);
		}
	}


    //substitution. What in detail?
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
	inline RetSeries BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSub(const PosTuple &posTuple, SubCaches &subCaches, const ArgsTuple &argsTuple) const
	{
        static_assert((boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value), "Positional and arguments' tuples' lengths do not match.");

		RetSeries retval;
		const const_iterator itf = end();
		for (const_iterator it = begin(); it != itf; ++it) 
		{
			RetSeries tmp = SubFunctor::template run<RetSeries>(it->cf, posTuple, subCaches, argsTuple);
			// NOTICE: series multadd here?
			tmp.baseMultBy(SubFunctor::template run<RetSeries>(it->key, posTuple, subCaches, argsTuple), argsTuple);

			retval.baseAdd(tmp, argsTuple);
		}

		return retval;
	}


    // empties the container of the BaseSeries object i.e. the container and the series is empty. is that a good idea????
    //it raises the question of 'what is an empty series?'
	template <PIRANHA_BASE_SERIES_TP_DECL>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::clearTerms()
	{
		container.clear();
	}


	// NOTE: can we use the concepts of next_echelon_type and echelon level here? Maybe we can avoid the runtime assert in numerical_cf?
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Series, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::baseSplit(std::vector<std::vector<Series> > &retval, const int n, const ArgsTuple &argsTuple) const
	{
		PIRANHA_ASSERT(retval.empty());
		PIRANHA_ASSERT(n >= 0 && n < boost::tuples::length<ArgsTuple>::value);

		if (n == 0) 
		{
			try {
				const std::vector<typename Derived::TermType const *> s(derived_const_cast->template get_sorted_series<Derived>(argsTuple));
				genericBaseSplit(retval, s.begin(), s.end(), argsTuple);

			} catch (const value_error &) 
			{
				genericBaseSplit(retval, begin(), end(), argsTuple);
			}

		} else 
		{
			if (!isSingleCf()) 
			{
				PIRANHA_THROW(value_error, "cannot split up to the specified level: series is non-degenerate");
			}

			begin()->cf.split(retval, n - 1, argsTuple);
		}
	}


	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class Series, class ArgsTuple>
	inline void BaseSeries<PIRANHA_BASE_SERIES_TP>::genericBaseSplit(std::vector<std::vector<Series> > &retval, const Iterator &start,
																		 const Iterator &end, const ArgsTuple &argsTuple) const
	{
		for (Iterator it = start; it != end; ++it) 
		{
			Series tmp_cf(Series::baseSeriesFromCf(FromIterator<Iterator>::get(it)->cf, argsTuple));
			Series tmp_key(Series::baseSeriesFromKey(FromIterator<Iterator>::get(it)->key, argsTuple));
			std::vector<Series> tmp;
			tmp.reserve(2);
			tmp.push_back(tmp_cf);
			tmp.push_back(tmp_key);

			retval.push_back(tmp);
		}
	}


	/// Return a vector of flattened terms.
	/**
	 * Flattened terms have coefficient series with a single term in all echelon levels.
	 */
	template <PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline std::vector<typename BaseSeries<PIRANHA_BASE_SERIES_TP>::TermType>
		BaseSeries<PIRANHA_BASE_SERIES_TP>::flattenTerms(const ArgsTuple &argsTuple) const
	{
		std::vector<TermType> retval;

		for (TermType term : *this)
		{
			SeriesFlattener<echelonLevel>::run(term.cf, term, retval, argsTuple);
		}
		return retval;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
