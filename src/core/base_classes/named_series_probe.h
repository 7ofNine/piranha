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

#ifndef PIRANHA_NAMED_SERIES_PROBE_H
#define PIRANHA_NAMED_SERIES_PROBE_H

#include <cstddef>
#include <utility>
#include <vector>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "../config.h" // For (un)likely
#include "../mp.h"
#include "../Psym.h"
#include "named_series_def.h"
#include "named_series_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
    //
	// TMP for arguments compatibility check.
    //
	template <class ArgsTuple>
	inline bool namedSeriesIsArgsCompatible(ArgsTuple const & argsTuple1, ArgsTuple const & argsTuple2)
	{
        // should we check that both ArgsTuples have the same number of components/vectors in them?
		const std::size_t w = argsTuple2.get_head().size(); // the size of the vector<Psym> at the head to the tuple of second parameter

		if (w > argsTuple1.get_head().size())               // not compatible if more arguments than in first tuple.
		{
			return false;
		}

		for (std::size_t i = 0; i < w; ++i) 
		{
			if (argsTuple1.get_head()[i] != argsTuple2.get_head()[i])  // same symbol in same position
			{
				return false;                                          // not the same symbol in same location they are not compatible
			}
		}

        // recurse into argsTuple(1/2). If we made it here we are compatible so far.
        // continue with the next vector in the ArgsTuples, i.e. tail
		return namedSeriesIsArgsCompatible(argsTuple1.get_tail(), argsTuple2.get_tail());
	}


    //
    // terminate recursion fur argument compatibility
    //
	template <>
	inline bool namedSeriesIsArgsCompatible<boost::tuples::null_type>(boost::tuples::null_type const &, boost::tuples::null_type const &)
	{
		return true;
	}


	/// Compatibility check for arguments.
	/**
	 * Test whether series' arguments are compatible with those from series2. Compatibility
	 * means that the number of arguments in all arguments sets/tuple elements are equal to or greater than series2's, and
	 * that arguments have the same positions as in series2's.
	 * @param[in] series2 series compatibility is tested against.
	 */
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline bool NamedSeries<PIRANHA_NAMED_SERIES_TP>::isArgsCompatible(Derived2 const &series2) const
	{
		// Use getter in second place because we may be interacting with other series type.
		return namedSeriesIsArgsCompatible(argumentsTuple, series2.arguments());
	}

    //
    // norm() of the series. norm is what?
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline double NamedSeries<PIRANHA_NAMED_SERIES_TP>::norm() const
	{
		return derived_const_cast->baseNorm(argumentsTuple);
	}


    //
    //  eval() : evaluate series at time t
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename TermEvalTypeDeterminer<Term>::Type
	NamedSeries<PIRANHA_NAMED_SERIES_TP>::eval(double const &t) const
	{
		return derived_const_cast->baseEval(t, argumentsTuple);
	}


    //
    // TMP function for checking that evaluation dictionary has all the elements needed.
    //
	template <class ArgsTuple>
	static inline bool checkEvalDict(EvalDict const &dictionary, ArgsTuple const &argsTuple)
	{
		std::size_t const size = argsTuple.get_head().size();
		EvalDict::const_iterator const itEnd = dictionary.end();

		for (std::size_t i = 0; i < size; ++i)
        {
			// If the dictionary does not contain the symbol's name, return false.
			if (dictionary.find(argsTuple.get_head()[i].getName()) == itEnd)
            {
				return false;
			}
		}

		// recurse to next tuple position.
		return checkEvalDict(dictionary, argsTuple.get_tail());
	}


    //
	// TMP function for checking that evaluation dictionary has all the elements needed.
	//
    // terminate recursion
    static inline bool checkEvalDict(EvalDict const &, boost::tuples::null_type const &)
	{
		return true;
	}


    //
    //
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename TermEvalTypeDeterminer<Term>::Type
	NamedSeries<PIRANHA_NAMED_SERIES_TP>::eval(EvalDict const &dictionary) const
	{
		if (!checkEvalDict(dictionary, argumentsTuple)) 
		{
			PIRANHA_THROW(value_error, "Evaluation dictionary does not contain entries for all the symbols of the series");
		}

		// Vector of original time evaluation vectors, paired with the corresponding symbol.
		// NOTE: here we allocate dynamically. This can be avoided by fixing a max number of items in time evaluation for psyms
		// and using static vectors. We should test performance before bothering though.
		std::vector<std::pair<Psym, std::vector<double> > > originalEval;
		originalEval.reserve(dictionary.size());

		const EvalDict::const_iterator itEnd(dictionary.end());

		for (EvalDict::const_iterator it = dictionary.begin(); it != itEnd; ++it)
        {
			originalEval.push_back(std::make_pair(Psym(it->first), Psym(it->first).getTimeEval()));
			Psym(it->first, it->second);
		}

		typename TermEvalTypeDeterminer<Term>::Type const retval(eval(0));

		// Restore original evaluation vectors.
		// NOTE: here RAII here, to be exception safe?
		for (std::size_t i = 0; i < originalEval.size(); ++i)
        {
			originalEval[i].first.setTimeEval(originalEval[i].second);
		}

		return retval;
	}


    //
    //  psi() determine power series iterations
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline int NamedSeries<PIRANHA_NAMED_SERIES_TP>::psi(int const start, int const step) const
	{
		return derived_const_cast->psIterations(start, step, argumentsTuple);
	}


    //
    // check if the tuple elements, which are vectors, are of the same size
    //
    // terminate recursion for the tuple Element size
	inline bool tupleVectorSameSizes(boost::tuples::null_type const &, boost::tuples::null_type const&)
	{
		return true;
	}


    //
    // check that the tuple elements are of same size (the elements are vectors)
    //
	template <class Tuple1, class Tuple2>
	inline bool tupleVectorSameSizes(const Tuple1 &tuple1, const Tuple2 &tuple2)
	{
		if (tuple1.get_head().size() != tuple2.get_head().size())  // here we test the size. get_head gives the element (i.e. vector)
        {
			return false;
		}

		return tupleVectorSameSizes(tuple1.get_tail(), tuple2.get_tail());
	}


    //
    // check on equality
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline bool NamedSeries<PIRANHA_NAMED_SERIES_TP>::isEqualTo(T const &other) const
	{
		// If the sizes of the arguments tuple elements do not coincide, series are different.
		if (!tupleVectorSameSizes(argumentsTuple, other.argumentsTuple)) 
        {
			return false;
		}

		// If arguments tuples are completely identical, just run the base
		// comparison function.
		if (argumentsTuple == other.arguments()) 
        {
			return derived_const_cast->baseEqualTo(other);
		}

        // different number of terms can not be equal
		if (derived_const_cast->length() != other.length())
        {
			return false;
		}

		// If we have same sizes of arguments tuples but they are not identical, then we may have to do
		// an arguments merge an see if the arguments are permutations of each other or if they are really different.
		// NOTE: this check is repeated in baseEqualTo, but doing it here could save a lot of work below.

		// Build a tuple of layouts.
		//typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type layout;
        LayoutTuple<typename Derived::ArgumentsDescription>::Type layoutTuple;
		// Get the relative layouts of this wrt other and put the result into layout.
		NamedSeriesGetLayout<ArgsTupleType>::run(argumentsTuple, other.arguments(), layoutTuple);

		// If the layout is bigger than the current ags tuple, it means that it is not a permutation,
		// there are different arguments in this and other. Hence we can return false.
		if (!tupleVectorSameSizes(argumentsTuple, layoutTuple)) 
        {
			return false;
		}

		// In this last case, the arguments are the same but they are ordered differently. Need to copy
		// this into new series with correct ordering and then do the comparison.
		// Build an empty retval and assign it the same arguments as this.
		Derived tmp;
		tmp.argumentsTuple = argumentsTuple;

		// Apply the layout to the arguments tuple of retval.
		NamedSeriesApplyLayoutToArgs<ArgsTupleType>::run(tmp.argumentsTuple, other.arguments(), layoutTuple);
		
        // Apply the layout to all terms of this and insert them into tmp.
		derived_const_cast->applyLayoutToTerms(layoutTuple, tmp.argumentsTuple, tmp);
		
        // Now we can perform the comparison between tmp and other.
		return tmp.baseEqualTo(other);
	}


    //
    //  compare ==
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline bool NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator==(T const &x) const
	{
		return NamedSeriesEqualitySelector<T>::run(*derived_const_cast, x);
	}


    //
    //  compare  !=
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline bool NamedSeries<PIRANHA_NAMED_SERIES_TP>::operator!=(T const &x) const
	{
		return !(operator==(x));
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
