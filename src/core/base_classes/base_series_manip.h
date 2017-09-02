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


	// TODO: update doc here.
	/// High-level insertion function.
	/**
	 * This function is used to insert terms into a series. It requires that the number of arguments
	 * of each element of the term is smaller or equal to the series',
	 * otherwise an assertion fails and the program aborts. base_pseries::mergeArgs,
	 * base_pseries::append_cf_args, base_pseries::append_trig_args, etc. can be used to add the needed arguments
	 * to the series.
	 *
	 * This function performs some checks and then calls llInsert.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool CanonicalCheck, bool Sign, class Term2, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term2, const ArgsTuple &argsTuple)
	{
        //convert the inserted term to the type of the term (TermType) of the series it is being inserted to
		TermConverter<Term2, TermType> convertedTerm(term2, argsTuple);

		// Make sure the appropriate routines for the management of arguments have been called.
		PIRANHA_ASSERT(convertedTerm.result.cf.isInsertable(argsTuple) && convertedTerm.result.key.isInsertable(argsTuple));

		TermType *newTerm(0); //NULL pointer
		if (convertedTerm.result.cf.needsPadding(argsTuple) || convertedTerm.result.key.needsPadding(argsTuple)) 
		{
			newTerm = TermType::allocator.allocate(1);
			TermType::allocator.construct(newTerm, convertedTerm.result); /// Shouldn't we just overload new???
			newTerm->cf.padRight(argsTuple);
			newTerm->key.padRight(argsTuple);
		}

		if (CanonicalCheck) 
		{
			if (!convertedTerm.result.isCanonical(argsTuple)) 
			{
				if (newTerm == 0) 
				{
					newTerm = TermType::allocator.allocate(1);
					TermType::allocator.construct(newTerm, convertedTerm.result);
				}
				newTerm->canonicalise(argsTuple);
			}
		}

		const TermType *insertTerm(0);
		if (newTerm) 
		{
			insertTerm = newTerm;

		} else 
		{
			insertTerm = &convertedTerm.result;
		}
		llInsert<Sign>(*insertTerm, argsTuple);
		
		if (newTerm) 
		{
			TermType::allocator.destroy(newTerm);   // should we use a helper for RAII?
			TermType::allocator.deallocate(newTerm, 1);
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Term2, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term, const ArgsTuple &argsTuple)
	{
		insert<true, true>(term, argsTuple);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::insertRange(const Iterator &begin, const Iterator &end, const ArgsTuple &argsTuple)
	{
		for (Iterator it = begin; it != end; ++it) 
		{
			insert(*it, argsTuple);
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<__PIRANHA_BASE_SERIES_TP>::findTerm(const TermType &term) const
	{
		return container.find(term);
	}


    //the low level insert routine. All the adjustments are out of the way. See the PIRANHA_ASSERT
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::llInsert(const TermType &term, const ArgsTuple &argsTuple)
	{
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

			// Check if the new coefficient can be ignored.
			if (it->cf.isIgnorable(argsTuple)) 
			{
				eraseTerm(it);
			}
		}
	}


	// Insert a new term into the series
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::termInsertNew(const TermType &term, const ArgsTuple &argsTuple)
	{
        //insert the term into the container
		std::pair<const_iterator, bool> res(container.insert(term));

		PIRANHA_ASSERT(res.second); // check that insert worked

        //if sign invert the coefficient
		if (!Sign) 
		{
			res.first->cf.invertSign(argsTuple);
		}
	}

    //remove term given by the iterator
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::eraseTerm(const const_iterator &it)
	{
		container.erase(it);
	}


	/// Swap the terms with another series.
	/**
	 * All terms get swapped.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSwap(Derived &series2)
	{
		PIRANHA_ASSERT(derived_cast != &series2);
		container.swap(series2.container);  // swap the series terms by swapping the terms container
	}

    //
	// Apply an arguments layout to all terms and insert them into retval.
    // for layout see named_series_manip.h
    // 
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class LayoutTuple, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::applyLayoutToTerms(LayoutTuple const &layoutTuple, ArgsTuple const &argsTuple, Derived &retval) const
	{
		const_iterator const itEnd = end();
		for (const_iterator it = begin(); it != itEnd; ++it) 
		{
			TermType term(*it);
			term.cf.applyLayout( layoutTuple, argsTuple);
			term.key.applyLayout(layoutTuple, argsTuple);

			retval.insert(term, argsTuple);
		}
	}


    //
    // test if we can trim some arguments from the argsTuple and terms
    // return the result as trimFlags(tuple of vector<bool>, false: trim i.e. the symbol in the argumentTuple is not used
    //
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::trimTestTerms(TrimFlags &trimFlags) const
	{
		const const_iterator itEnd = end();
		for (const_iterator it = begin(); it != itEnd; ++it) 
		{
			it->cf.trimTest(trimFlags);
			it->key.trimTest(trimFlags);
		}
	}


    //
    // trim the arguments according to the trimFlags from the *this terms argumentTuple and return the modified series in retval
    // the triiming of the argumentsTuple was done earlier
    //
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::trimTerms(TrimFlags const &trimFlags, ArgsTuple const &argsTuple, Derived &retval) const
	{
		const const_iterator itEnd = end();
		for (const_iterator it = begin(); it != itEnd; ++it) 
		{
			retval.insert(TermType(typename TermType::CfType(it->cf.trim(trimFlags, argsTuple)), typename TermType::KeyType(it->key.trim(trimFlags, argsTuple))), argsTuple);
		}
	}


    //substitution. What in detail?
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
	inline RetSeries BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSub(const PosTuple &posTuple, SubCaches &subCaches, const ArgsTuple &argsTuple) const
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
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::clearTerms()
	{
		container.clear();
	}


	// NOTE: can we use the concepts of next_echelon_type and echelon level here? Maybe we can avoid the runtime assert in numerical_cf?
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Series, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSplit(std::vector<std::vector<Series> > &retval, const int n, const ArgsTuple &argsTuple) const
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


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Iterator, class Series, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::genericBaseSplit(std::vector<std::vector<Series> > &retval, const Iterator &start,
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
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class ArgsTuple>
	inline std::vector<typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::TermType>
		BaseSeries<__PIRANHA_BASE_SERIES_TP>::flattenTerms(const ArgsTuple &argsTuple) const
	{
		std::vector<TermType> retval;
		TermType term;
		const const_iterator itf(end());
		for (const_iterator it = begin(); it != itf; ++it) 
		{
			// Create the term that will be inserted at the end of the recursion.
			term = *it;
			// Start the recursion.
			SeriesFlattener<echelonLevel>::run(term.cf, term, retval, argsTuple);
		}
		return retval;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
