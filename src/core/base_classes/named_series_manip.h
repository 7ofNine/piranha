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

#ifndef PIRANHA_NAMED_SERIES_MANIP_H
#define PIRANHA_NAMED_SERIES_MANIP_H

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <boost/tuple/tuple.hpp> 

#include "../config.h" // For (un)likely
#include "../exceptions.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "named_series_def.h"

#define derivedConstCast static_cast<Derived const *>(this)
#define derivedCast      static_cast<Derived       *>(this)

namespace piranha
{
	// Template metaprogramming for swapping of arguments in named series.
	template <class ArgsTuple>
	class NamedSeriesSwap
    {
	    public: 
        	
        static void run(ArgsTuple &argsTuple1, ArgsTuple &argsTuple2)
        {
			argsTuple1.get_head().swap(argsTuple2.get_head());   // swap the psym vector tuple elements

			NamedSeriesSwap<typename ArgsTuple::tail_type>::run(argsTuple1.get_tail(), argsTuple2.get_tail()); //recurse to the next tuple element
		}
	};

    //terminate  recursion through tuple elements
	template <>
	class NamedSeriesSwap<boost::tuples::null_type>
    {
        public: 

		static void run(boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


	/// Swap contents of series. Avoids creation (performance!)
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::swap(Derived &series)
	{
		// Do something only if we are not swapping with self.
		if (derivedCast != &series)
        {
			NamedSeriesSwap<ArgsTupleType>::run(argumentsTuple, series.argumentsTuple); // swap the argument tuples (named series only)
			derivedCast->baseSwap(series); // swap the terms
		}
	}


	// Meta-programming for appending an argument.
    //
    // append the symbol "arg" to the proper "argsTuple" element  given by the "argumentType" (poly, trig)
    //
	template <class ArgsDescr>
	class NamedSeriesAppendArg
    {
        
        typedef typename NTuple<VectorPsym, boost::tuples::length<ArgsDescr>::value>::Type ArgsTupleType;

        public:

		static void run(std::string const &argumentType, ArgsTupleType &argsTuple, Psym const &arg)
		{
			if (ArgsDescr::head_type::name == argumentType) 
			{
				// Check that the argument is not already present in this set.
				for (VectorPsym::iterator it = argsTuple.get_head().begin(); it != argsTuple.get_head().end(); ++it) 
				{
					if (arg == (*it)) 
					{
						std::cout << "Error: " << argumentType << " argument '" << it->getName() << "' already present in the set.\n";
						return;
					}
				}
				argsTuple.get_head().push_back(arg);  // add new symbol to the end (it is a vector)

			} else 
			{
                // try next tuple element
				NamedSeriesAppendArg<typename ArgsDescr::tail_type>::run(argumentType, argsTuple.get_tail(), arg);
			}
		}
	};

    //terminate recursion for append argument
	template <>
	class NamedSeriesAppendArg<boost::tuples::null_type>
    {

        public: 

		static void run(std::string const &argumentType, boost::tuples::null_type const &, Psym const &) 
		{
			std::cout << "Error: '" << argumentType << "' arguments are not known." << std::endl;
		}
	};


    //
    // append argument/symbol by argument type (as defined in the argument descriptors) ("poly", "trig", ...)
	//
    // only to be done for an empty series i.e. during construction.
    // This also means that arguments are added to the series before the terms are added
    //
    template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::appendArg(std::string const &argumentType, Psym const &arg)
	{
		PIRANHA_ASSERT(derivedConstCast->empty());  // check empty series

		NamedSeriesAppendArg<ArgumentsDescription>::run(argumentType, argumentsTuple, arg);
	}

    //
    //apend argument/symbol at ecehelon level N in the argumentsTuple. See Arguments descriptor
    //
    //only to be done for an empty series i.e. during construction.
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::appendArg(Psym const &arg)
	{
        static_assert(N >= 0, "Trying to append argument with invalid index.");
		PIRANHA_ASSERT(derivedConstCast->empty());  // check empty series

		// Check that the argument is not already present in this set.
		for (VectorPsym::iterator it = argumentsTuple.template get<N>().begin(); it != argumentsTuple.template get<N>().end(); ++it) 
		{
			if (arg == (*it)) 
			{
				std::cout << "Error: argument '" << it->getName() << "' already present in the set." << std::endl;
				return;
			}
		}

		argumentsTuple.template get<N>().push_back(arg); // insert new symbol into ntuple element at position N.
	}


    //
    // combine two sets of arguments into one set
    //
    // the result is a modified argsTuple in *this that is the a sequenced union from the
    // *this and series2. The sequence of the resulting arguments is different from the sequences of the incoming arguments.
    // Also the result is dependent on which series is first, i.e. the result is non commutative
    // Merging with oneself is not allowed on this level. should we fail it?? This would be programming error
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::mergeArgs(Derived2 const &series2)
	{
		if (static_cast<void *>(this) == static_cast<void const *>(&series2)) 
		{
			PIRANHA_DEBUG(std::cout << "Trying to merge with self, returning." << std::endl);
			return;
		}

        // the same symbols have to be in the same location, if not they are not compatible
		if (!isArgsCompatible(series2))  // if they are not compatible i.e. argsTuple2 symbols are not all contained in argsTuple of this series. 
		{
			mergeIncompatibleArgs(series2); // create union of the symbols. Beware the sequence of symbols gets changed, not just appended.
		}
	}


	// Template metaprogramming for getting the relative layout of two series.
	// used during multiplication. 
    //
    //
    // Layout: are tuples with respect to the exponential and trigonometric keys as they are described by the argsTuple members of the series.
    // The elements of the tuple are vectors of pairs of (bool, int).
    // The layout is determined by starting in series2. If the symbol is present in series2 but not in series1 the
    // bool flag is set to false. The index int is not used (typically 0). The index of the pair in the vector determines which index into the 
    // argumentsTuple Psym vector it corresponds to (this could be changed to make it homogenous, see below)
    // If the symbol is present in series1 the bool is set to true and the index integer is set to the index in the argsTuple psym vector of series1.
    // The same is done for symbols that are present in series1 but not in series2.
    // e.g.
    // series2          :  ("x", "y", "z")
    // series1 (= *this):  ("u", "x", "z")
    //               "x"        "y"         "z"        "u"
    // => layout ((true, 1), (false, 0), (true, 2), (true,0))
    //
	template <class ArgsTuple>
	class NamedSeriesGetLayout
    {
       public:
       // typedef std::pair<bool, std::size_t> LayoutElement;
       // typedef std::vector<LayoutElement> Layout;
       // typedef typename NTuple< Layout, boost::tuples::length<ArgsTuple>::value >::Type LayoutTuple;

		static void run(ArgsTuple const &argsTuple1, ArgsTuple const &argsTuple2, typename LayoutTuple<ArgsTuple>::Type &layoutTuple) 
		{
			// Store frequently-used variables.
			VectorPsym const &symbols1   = argsTuple1.get_head(); 
			VectorPsym const &symbols2   = argsTuple2.get_head();
			std::size_t const size1      = symbols1.size(); 
			std::size_t const size2      = symbols2.size();
			Layout &layout               = layoutTuple.get_head();

			// First we must construct symbols2's layout wrt to symbols1.
			layout.resize(size2);
			for (std::size_t i = 0; i < size2; ++i) 
            {
				// Find current symbols2's symbol within symbols1.
				const VectorPsym::const_iterator result = std::find(symbols1.begin(), symbols1.end(), symbols2[i]);
				if (result == symbols1.end())  // not found
                {
					layout[i].first = false;
				} else 
                {
					// Found, mark and store its position in symbols1
					layout[i].first  = true;
					layout[i].second = result - symbols1.begin();
				}
			}

			// Handle those elements of symbols1 that are not yet present in
			// the layout (i.e., they are not in symbols2)
			for (std::size_t i = 0; i < size1; ++i) 
            {
				// Look for element index i in the layout.
				bool found = false;
				std::size_t const lSize = layout.size();
				for (std::size_t j = 0; j < lSize; ++j) 
                {
					if (layout[j].first && layout[j].second == i) 
                    {
                        //already handled by loop over symbols2
						found = true;
						break;
					}
				}

				// If we did not find it, append it to the layout.
                // It is in symbols1 but not in symbols2
				if (!found) 
                {
					layout.push_back(LayoutElement(true, i));
				}
			}

            // recurse to next argsTuple element
			NamedSeriesGetLayout<typename ArgsTuple::tail_type>::run(argsTuple1.get_tail(), argsTuple2.get_tail(), layoutTuple.get_tail());
		}
	};


    //
    // terminate recursion for determining layout
    //
	template <>
	class NamedSeriesGetLayout<boost::tuples::null_type>
	{
		public:

		static void run(boost::tuples::null_type const &, boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


    //
	// Template metaprogramming for applying a layout to a series.
    // modify arguments for argsTuple1 with argsTuple2 according to the layout
    // TODo: example
    //
	template <class ArgsTuple>
	class NamedSeriesApplyLayoutToArgs
    {
        public:

		static void run(ArgsTuple &argsTuple1, ArgsTuple const &argsTuple2, typename LayoutTuple<ArgsTuple>::Type const &layoutTuple) 
		{
			// Store frequently-used variables.
			VectorPsym        &symbols1 = argsTuple1.get_head();
			VectorPsym const  &symbols2 = argsTuple2.get_head();
			Layout const      &layout   = layoutTuple.get_head();
			std::size_t const  lSize    = layout.size();
			
            // The layout must have at least all arguments in symbols1.
			PIRANHA_ASSERT(lSize >= symbols1.size());
			
            // Memorize the old vector.
			const VectorPsym oldSymbols(symbols1);
			// Make space.
			symbols1.reserve(lSize);

			for (std::size_t i = 0; i < lSize; ++i) 
            {
				if (layout[i].first) 
                {
					// The argument was present in the old arguments sets. Copy it over.
					PIRANHA_ASSERT(layout[i].second < oldSymbols.size());
					if (i < symbols1.size())
					{
						symbols1[i] = oldSymbols[layout[i].second];

					} else 
					{
						symbols1.push_back(oldSymbols[layout[i].second]);
					}
				} else 
                {
					// The argument was not present in the old arguments set (symbols1). Fetch it from symbols2.
					PIRANHA_ASSERT(i < symbols2.size());
					if (i < symbols1.size()) 
                    {
						symbols1[i] = symbols2[i];

					} else 
                    {
						symbols1.push_back(symbols2[i]);
					}
				}
			}

            // recurse into the next argsTuple element
			NamedSeriesApplyLayoutToArgs<typename ArgsTuple::tail_type>::run(argsTuple1.get_tail(), argsTuple2.get_tail(), layoutTuple.get_tail());
		}
	};


    //
    // terminate recursion for applying layout to argsTuple
	template <>
	class NamedSeriesApplyLayoutToArgs<boost::tuples::null_type>
	{
		public:

			static void run(boost::tuples::null_type const &, boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


    //
    // Layout: are tuples with respect to the exponential and trigonometric keys as they are described by the argsTuple members of the series.
    // The elements of the layout tuple are vectors of pairs of (bool, int).
    // The layout is determined by starting in series2. If the symbol is present in series2 but not in series1 the
    // bool flag is set to false. The index int is not used (typically 0). The index of the pair in the vector determines which index into the 
    // argumentsTuple Psym vector of series 2 it corresponds to (this could be changed to make it homogenous, see below)
    // If the symbol is present in series1 the bool is set to true and the index integer is set to the index in the argsTuple psym vector of series1.
    // The same is done for symbols that are present in series1 but not in series2.
    // e.g.
    // series2          :  ("x", "y", "z")
    // series1 (= *this):  ("u", "x", "z")
    //               "x"        "y"         "z"        "u"
    // => layout ((true, 1), (false, 0), (true, 2), (true,0))
    //
    // here we create the union of the two argTuples and update *this acordingly
    //
    // note that the merged arguments start with the arguments in sequence as they derive from series 2(!) followed in sequence by the still missing arguments from series 1.
    // This is explicitly used in the code for multiplication/division etc. I.e. within the code, series 1 and series 2 (and its terms ) don't commute!. This way terms of series 2 can 
    // be directly merged into the resulting series, their arguments are the starting sequence of the resulting arguments sequence.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::mergeIncompatibleArgs(Derived2 const &series2)
	{
		// Build an empty retval and assign to it the same arguments as this.
		Derived retval;
		retval.argumentsTuple = argumentsTuple;

		// Build layout tuple.
		//typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type layout;
        typename LayoutTuple<ArgsTupleType>::Type layoutTuple;
		// Get the relative layouts of this wrt series2 and put the result into layout.
		NamedSeriesGetLayout<ArgsTupleType>::run(retval.argumentsTuple, series2.arguments(), layoutTuple);
		
		// Apply the layout to the arguments tuple of retval.
		NamedSeriesApplyLayoutToArgs<ArgsTupleType>::run(retval.argumentsTuple, series2.arguments(), layoutTuple);
		
		// Apply the layout to all terms of this, which will be inserted into retval.
		derivedConstCast->applyLayoutToTerms(layoutTuple, retval.argumentsTuple, retval);
		
		// Finally, swap the contents of retval with this. this is now a series with argumentTuple that combines this and series2
        // i.e. operations like addition and multiplication can now be executed.
		swap(retval);
	}

    //
    // set up flags tuple for triming unneeded arguments
    // tuple elements are vectors<bool> 
    // false: remove  (argument is not present)
    // true:  don't remove (argument is present)
    //
	// this is just a complicated way to set all the needed flags in there appropriate vectors in the tuples to false

	template <class TrimFlags, class ArgsTuple>
	class TrimFlagsInit
    {
        public:

		static void run(TrimFlags &trimFlags, ArgsTuple const &argsTuple)
        {
			std::size_t const size = argsTuple.get_head().size();
			trimFlags.get_head().resize(size);
			for (std::size_t i = 0; i < size; ++i) 
            {
				trimFlags.get_head()[i] = false;
			}

			TrimFlagsInit<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(trimFlags.get_tail(), argsTuple.get_tail());
		}
	};


    //
    // terminate trim flag tuple initialization
    //
	template <>
	class TrimFlagsInit<boost::tuples::null_type, boost::tuples::null_type>
    {
        public:

		static void run(boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


	inline bool trimFlagsProceed(boost::tuples::null_type const &)
	{
		return false;
	}

    //
    // check if there is something to trim from the arguments
    // TrimFlags is a tuple of vector<bool> elements
    //
	template <class TrimFlags>
	inline bool trimFlagsProceed(TrimFlags const &trimFlags)
	{
		std::size_t const size = trimFlags.get_head().size();

		for (std::size_t i = 0; i < size; ++i) 
		{
			// If we find a flag that was never turned on, we have something to trim.
			if (!trimFlags.get_head()[i]) 
			{
				return true;
			}
		}

        // recurse to next tuple
		return trimFlagsProceed(trimFlags.get_tail());
	}


    //
    // trim argumentTuple by the arguments indicated by the trimflagsTuple
    //
	template <class TrimFlags, class ArgsTuple>
	class TrimArguments 
	{
        public:

		static void run(TrimFlags const &trimFlags, ArgsTuple &argsTuple) 
		{
			const std::size_t size = trimFlags.get_head().size();

			PIRANHA_ASSERT(size == argsTuple.get_head().size());
			
            VectorPsym newSymbols;
			for (std::size_t i = 0; i < size; ++i) 
			{
				if (trimFlags.get_head()[i]) // if false we do't push i.e. argument is removed.
				{
					newSymbols.push_back(argsTuple.get_head()[i]);
				}
			}

			newSymbols.swap(argsTuple.get_head()); // replace the old arguments with the new arguments (vector.swap()).

            // recurse to the next argumentsTuple element, i.e. key in the echelon level
			TrimArguments<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(trimFlags.get_tail(), argsTuple.get_tail());
		}
	};


    //
    // terminate recursion for argument trimming at end of argumentsTuple
    //
	template <>
	class TrimArguments<boost::tuples::null_type, boost::tuples::null_type> 
	{
        public:

		static void run(boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


    //
    // remove arguments from argumentTuple that are no longer needed (i.e. evaluate to zero or 1)
    // and adjust *this accordingly
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::trim()
	{
		typedef typename NTuple<std::vector<bool>, Derived::echelonLevel + 1>::Type TrimFlagsType; 
		TrimFlagsType trimFlags;

		TrimFlagsInit<TrimFlagsType, ArgsTupleType>::run(trimFlags, argumentsTuple);
		
        derivedConstCast->trimTestTerms(trimFlags); // check if arguments can be trimmed

		if (trimFlagsProceed(trimFlags))   // something to trim
		{
			// First let's do the arguments.
			TrimArguments<TrimFlagsType, ArgsTupleType>::run(trimFlags, argumentsTuple); // this changes the arguments tuple of this series
			// Let's proceed to the terms now.                                           // this series goes temporarily into an inconsistent   
			Derived tmpSeries;                                                           // state until baseSwap
			derivedCast->trimTerms(trimFlags, argumentsTuple, tmpSeries);                // trimTerms is in BaseSeries 
			derivedCast->baseSwap(tmpSeries);   // exchange trimmed contents (== tmpSeries) with *this. i.e. update, argumentsTuple is already done
		}
	}

    
    //TODO: how does substitution work??
    //
    // Initialization functor for substitution caches
	//
    template <class SubstitutionCaches, class SubstitutionSeries, class ArgsTuple>
	class InitSubstitutionCaches
	{
        public:

		static void run(SubstitutionCaches &substitutionCaches, SubstitutionSeries const &series, ArgsTuple const *argsTuple) 
		{
			substitutionCaches.get_head().setup(series, argsTuple);
			InitSubstitutionCaches<typename SubstitutionCaches::tail_type, SubstitutionSeries, ArgsTuple>::run(substitutionCaches.get_tail(), series, argsTuple);
		}
	};


    //
    // terminate recursion for substitution caches intitializtion
    //
	template <class SubstitutionSeries, class ArgsTuple>
	class InitSubstitutionCaches<boost::tuples::null_type, SubstitutionSeries, ArgsTuple>
	{
        public:

		static void run(boost::tuples::null_type const &, SubstitutionSeries const &, ArgsTuple const *) {}
	};


    //
    // substitute series for argument with name 
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class SubstitutionSeries>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::sub(std::string const &name, SubstitutionSeries const &series) const
	{
     // How to test that "name" is actually present in ArgsTuple and we can substitute it?? We already had that problem that substitution didn't fail
        typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type    PositionTuple;
        PositionTuple const testPositionTuple = psyms2pos(VectorPsym(1, Psym(name)), argumentsTuple);
        
		typedef typename Derived::TermType::CfType::
			template SubstitutionCacheSelector<SubstitutionSeries, typename Derived::TermType::KeyType::
			template SubstitutionCacheSelector<SubstitutionSeries, boost::tuples::null_type, ArgsTupleType>::Type, ArgsTupleType>::Type    SubstitutionCaches;

		typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type    PositionTuple;

        static_assert(boost::tuples::length<SubstitutionCaches>::value == boost::tuples::length<PositionTuple>::value,
			"Size mismatch for position and cache tuples in series substitution.");

		SubstitutionCaches substitutionCaches;

		Derived original(*derivedConstCast);

		SubstitutionSeries substitutionCopy(series);
		original.mergeArgs(substitutionCopy); // we already merge arguments but we don't know yet if argument with "name" actually exists.
        substitutionCopy.mergeArgs(original);

		// Init substitution caches using s_copy and this_copy.m_arguments.
		InitSubstitutionCaches<SubstitutionCaches, SubstitutionSeries, ArgsTupleType>::run(substitutionCaches, substitutionCopy, &original.argumentsTuple);

		PositionTuple const positionTuple = psyms2pos(VectorPsym(1, Psym(name)), original.argumentsTuple);

		Derived retval(original.template baseSub<Derived, typename Derived::SubstitutionFunctor>(positionTuple, substitutionCaches, original.argumentsTuple));

		retval.argumentsTuple = original.argumentsTuple;
		retval.trim();

		return retval;
	}

    //
    // Split does what?
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<std::vector<Derived> > NamedSeries<PIRANHA_NAMED_SERIES_TP>::split(const int n) const
	{
		if (n < 0 || n >= boost::tuples::length<ArgsTupleType>::value) 
		{
			PIRANHA_THROW(value_error, "Splitting level must be a non-negative integer less than the echelon level of the series");
		}

		std::vector<std::vector<Derived> > retval;
		derivedConstCast->baseSplit(retval, n, argumentsTuple);

		std::size_t const size = retval.size();
		for (std::size_t i = 0; i < size; ++i) 
		{
			retval[i][0].argumentsTuple = argumentsTuple;
			retval[i][0].trim();

			retval[i][1].argumentsTuple = argumentsTuple;
			retval[i][1].trim();
		}

		return retval;
	}


	// Flatten series.
	/**
	 * At echelon level 0, this method will return a vector of series each one consisting of a single term of the original series.
	 * For echelon levels above 0, each coefficient series will also consist of only one term. This method, in other words, flattens
	 * out the tree structure of echeloned series into a vector of single-term series, whose coefficients series are also single-term series
	 * and so on. For instance, the Poisson series
	 * \f[
	 * \left(x + y\right) \cos a + z
	 * \f]
	 * will be flattened into the vector \f$ \left[ x\cos a, y\cos a, z \right] \f$.
	 */
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<Derived> NamedSeries<PIRANHA_NAMED_SERIES_TP>::flatten() const
	{
		std::vector<typename Derived::TermType> const tmp(derivedConstCast->flattenTerms(argumentsTuple));
		std::vector<Derived> retval;

        typedef typename std::vector<typename Derived::TermType>::const_iterator Iterator;
		Iterator const itEnd(tmp.end());

		for  (Iterator it = tmp.begin(); it != itEnd; ++it) 
		{
			Derived tmpSeries;
			tmpSeries.insert(*it, argumentsTuple);
			tmpSeries.argumentsTuple = argumentsTuple;
			tmpSeries.trim();

			retval.push_back(Derived());
			retval.back().swap(tmpSeries);
		}

		return retval;
	}
}

#undef derivedConstCast
#undef derivedCast

#endif
