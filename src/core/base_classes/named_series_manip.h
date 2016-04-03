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
			argsTuple1.get_head().swap(argsTuple2.get_head());

			NamedSeriesSwap<typename ArgsTuple::tail_type>::run(argsTuple1.get_tail(), argsTuple2.get_tail());
		}
	};


	template <>
	class NamedSeriesSwap<boost::tuples::null_type>
    {
        public: 

		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	/// Swap contents of series.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::swap(Derived &series)
	{
		// Do something only if we are not swapping with self.
		if (derivedCast != &series)
        {
			NamedSeriesSwap<ArgsTupleType>::run(argumentsTuple, series.argumentsTuple);
			derivedCast->baseSwap(series);
		}
	}


	// Meta-programming for appending an argument.
	template <class ArgsDescr>
	class NamedSeriesAppendArg
    {

        public:

		static void run(std::string const &s, typename NTuple<VectorPsym, boost::tuples::length<ArgsDescr>::value>::Type &argsTuple, Psym const &arg)
		{
			if (ArgsDescr::head_type::name == s) 
			{
				// Check that the argument is not already present in this set.
				for (VectorPsym::iterator it = argsTuple.get_head().begin(); it != argsTuple.get_head().end(); ++it) 
				{
					if (arg == (*it)) 
					{
						std::cout << "Error: " << s << " argument '" << it->getName() << "' already present in the set.\n";
						return;
					}
				}
				argsTuple.get_head().push_back(arg);

			} else 
			{
				NamedSeriesAppendArg<typename ArgsDescr::tail_type>::run(s, argsTuple.get_tail(), arg);
			}
		}
	};


	template <>
	class NamedSeriesAppendArg<boost::tuples::null_type>
    {

        public: 

		static void run(std::string const &s, boost::tuples::null_type const &, Psym const &) 
		{
			std::cout << "Error: '" << s << "' arguments are not known." << std::endl;
		}
	};


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::appendArg(std::string const &s, Psym const &arg)
	{
		PIRANHA_ASSERT(derivedConstCast->empty());

		NamedSeriesAppendArg<ArgumentsDescription>::run(s, argumentsTuple, arg);
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::appendArg(Psym const &arg)
	{
		PIRANHA_STATIC_CHECK(N >= 0, "Trying to append argument with invalid index.");
		PIRANHA_ASSERT(derivedConstCast->empty());

		// Check that the argument is not already present in this set.
		for (VectorPsym::iterator it = argumentsTuple.template get<N>().begin(); it != argumentsTuple.template get<N>().end(); ++it) 
		{
			if (arg == (*it)) 
			{
				std::cout << "Error: argument '" << it->getName() << "' already present in the set." << std::endl;
				return;
			}
		}

		argumentsTuple.template get<N>().push_back(arg);
	}



	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::mergeArgs(Derived2 const &series2)
	{
		if (static_cast<void *>(this) == static_cast<void const *>(&series2)) 
		{
			PIRANHA_DEBUG(std::cout << "Trying to merge with self, returning." << std::endl);
			return;
		}

		if (unlikely(!isArgsCompatible(series2))) 
		{
			mergeIncompatibleArgs(series2);
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
        typedef std::pair<bool, std::size_t> LayoutElement;
        typedef std::vector<LayoutElement> Layout;
        typedef typename NTuple< Layout, boost::tuples::length<ArgsTuple>::value >::Type LayoutTuple;

		static void run(ArgsTuple const &argsTuple1, ArgsTuple const &argsTuple2, LayoutTuple &layoutTuple) 
		{
			// Store frequently-used variables.
			VectorPsym const &symbols1   = argsTuple1.get_head(); 
			VectorPsym const &symbols2   = argsTuple2.get_head();
			std::size_t const size1      = symbols1.size(); 
			std::size_t const size2      = symbols2.size();
			Layout &layout = layoutTuple.get_head();
			// First we must construct symbols2's layout wrt to symbols1.
			layout.resize(size2);
			for (std::size_t i = 0; i < size2; ++i) 
            {
				// Let's see if current symbols2's symbol is present in symbols1.
				const VectorPsym::const_iterator result = std::find(symbols1.begin(), symbols1.end(), symbols2[i]);
				if (result == symbols1.end()) 
                {
					layout[i].first = false;
				} else 
                {
					// If present, mark its position in symbols1
					layout[i].first  = true;
					layout[i].second = result - symbols1.begin();
				}
			}

			// Now we must take care of those elements of symbols1 that are not represented in
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

			NamedSeriesGetLayout<typename ArgsTuple::tail_type>::run(argsTuple1.get_tail(), argsTuple2.get_tail(), layoutTuple.get_tail());
		}
	};


	template <>
	class NamedSeriesGetLayout<boost::tuples::null_type>
	{
		public:

		static void run(boost::tuples::null_type const &, boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};



	// Template metaprogramming for applying a layout to a series.
	template <class ArgsTuple>
	class NamedSeriesApplyLayoutToArgs
    {
        public:

		static void run(ArgsTuple &argsTuple1, ArgsTuple const &argsTuple2, typename NTuple < std::vector<std::pair<bool, std::size_t> >,
						boost::tuples::length<ArgsTuple>::value >::Type const &layout) 
		{
			// Store frequently-used variables.
			VectorPsym       &symbols1 = argsTuple1.get_head();
			const VectorPsym &symbols2 = argsTuple2.get_head();
			const std::vector<std::pair<bool, std::size_t> > &l = layout.get_head();
			const std::size_t l_size = l.size();
			// The layout must have at least all arguments in v1.
			PIRANHA_ASSERT(l_size >= symbols1.size());
			// Memorize the old vector.
			const VectorPsym old(symbols1);
			// Make space.
			symbols1.reserve(l_size);

			for (std::size_t i = 0; i < l_size; ++i) 
            {
				if (l[i].first) 
                {
					// The argument was present in the old arguments sets. Copy it over.
					PIRANHA_ASSERT(l[i].second < old.size());
					if (i < symbols1.size())
					{
						symbols1[i] = old[l[i].second];

					} else 
					{
						symbols1.push_back(old[l[i].second]);
					}
				} else 
                {
					// The argument was not present in the old arguments sets. Fetch it from a2.
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

			NamedSeriesApplyLayoutToArgs<typename ArgsTuple::tail_type>::run(argsTuple1.get_tail(), argsTuple2.get_tail(), layout.get_tail());
		}
	};


	template <>
	class NamedSeriesApplyLayoutToArgs<boost::tuples::null_type>
	{
		public:

			static void run(boost::tuples::null_type const &, boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


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
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::mergeIncompatibleArgs(Derived2 const &series2)
	{
		// Build an empty retval and assign it the same arguments as this.
		Derived retval;
		retval.argumentsTuple = argumentsTuple;

		// Build layout tuple.
		//typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type layout;
        typename NamedSeriesGetLayout<ArgsTupleType>::LayoutTuple layoutTuple;
		// Get the relative layouts of this wrt series2 and put the result into layout.
		NamedSeriesGetLayout<ArgsTupleType>::run(retval.argumentsTuple, series2.arguments(), layoutTuple);
		
		// Apply the layout to the arguments tuple of retval.
		NamedSeriesApplyLayoutToArgs<ArgsTupleType>::run(retval.argumentsTuple, series2.arguments(), layoutTuple);
		
		// Apply the layout to all terms of this, which will be inserted into retval.
		derivedConstCast->applyLayoutToTerms(layoutTuple, retval, retval.argumentsTuple);
		
		// Finally, swap the contents of retval with this.
		swap(retval);
	}


	template <class TrimFlags, class ArgsTuple>
	class TrimFlagsInit
    {
        public:

		static void run(TrimFlags &tf, const ArgsTuple &argsTuple)
        {
			const std::size_t size = argsTuple.get_head().size();
			tf.get_head().resize(size);
			for (std::size_t i = 0; i < size; ++i) 
            {
				tf.get_head()[i] = false;
			}

			TrimFlagsInit<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(tf.get_tail(), argsTuple.get_tail());
		}
	};


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


	template <class TrimFlags>
	inline bool trimFlagsProceed(TrimFlags const &tf)
	{
		std::size_t const size = tf.get_head().size();

		for (std::size_t i = 0; i < size; ++i) 
		{
			// If we find a flag that was never turned on, we have something to trim.
			if (!tf.get_head()[i]) 
			{
				return true;
			}
		}

		return trimFlagsProceed(tf.get_tail());
	}


	template <class TrimFlags, class ArgsTuple>
	class TrimArguments 
	{
        public:

		static void run(const TrimFlags &tf, ArgsTuple &argsTuple) 
		{
			const std::size_t size = tf.get_head().size();
			PIRANHA_ASSERT(size == argsTuple.get_head().size());
			VectorPsym new_vector;
			for (std::size_t i = 0; i < size; ++i) 
			{
				if (tf.get_head()[i]) 
				{
					new_vector.push_back(argsTuple.get_head()[i]);
				}
			}

			new_vector.swap(argsTuple.get_head());
			TrimArguments<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(tf.get_tail(), argsTuple.get_tail());
		}
	};


	template <>
	class TrimArguments<boost::tuples::null_type, boost::tuples::null_type> 
	{
        public:

		static void run(boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::trim()
	{
		typedef typename NTuple<std::vector<char>, Derived::echelonLevel + 1>::Type TrimFlagsType;
		TrimFlagsType trimFlags;

		TrimFlagsInit<TrimFlagsType, ArgsTupleType>::run(trimFlags, argumentsTuple);
		
        derivedConstCast->trimTestTerms(trimFlags);

		if (trimFlagsProceed(trimFlags)) 
		{
			// First let's do the arguments.
			TrimArguments<TrimFlagsType, ArgsTupleType>::run(trimFlags, argumentsTuple);
			// Let's proceed to the terms now.
			Derived tmpSeries;
			derivedCast->trimTerms(trimFlags, tmpSeries, argumentsTuple);
			derivedCast->baseSwap(tmpSeries);
		}
	}

    // Initialization functor for substitution cache
	template <class SubCaches, class SubSeries, class ArgsTuple>
	class InitSubCaches
	{
        public:

		static void run(SubCaches &sub_caches, SubSeries const &s, ArgsTuple const *argsTuple) 
		{
			sub_caches.get_head().setup(s, argsTuple);
			InitSubCaches<typename SubCaches::tail_type, SubSeries,ArgsTuple>::run(sub_caches.get_tail(), s, argsTuple);
		}
	};


	template <class SubSeries, class ArgsTuple>
	class InitSubCaches<boost::tuples::null_type, SubSeries, ArgsTuple>
	{
        public:

		static void run(boost::tuples::null_type const &, SubSeries const &, ArgsTuple const *) {}
	};

    //
    // substitute series s for argument name name
    //
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class SubSeries>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::sub(const std::string &name, const SubSeries &s) const
	{
		typedef typename Derived::TermType::CfType::
			template SubstitutionCacheSelector<SubSeries, typename Derived::TermType::KeyType::
			template SubstitutionCacheSelector<SubSeries, boost::tuples::null_type, ArgsTupleType>
			::type, ArgsTupleType>::type    sub_caches_type;

		typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type    pos_tuple_type;

		PIRANHA_STATIC_CHECK(boost::tuples::length<sub_caches_type>::value == boost::tuples::length<pos_tuple_type>::value,
			"Size mismatch for position and cache tuples in series substitution.");

		const Psym p(name);
		sub_caches_type sub_caches;
		Derived this_copy(*derivedConstCast);
		SubSeries s_copy(s);
		this_copy.mergeArgs(s_copy);
		s_copy.mergeArgs(this_copy);

		// Init substitution caches using s_copy and this_copy.m_arguments.
		InitSubCaches<sub_caches_type, SubSeries, ArgsTupleType>::run(sub_caches, s_copy, &this_copy.argumentsTuple);

		const pos_tuple_type pos_tuple = psyms2pos(VectorPsym(1, p), this_copy.argumentsTuple);

		Derived retval(this_copy.template baseSub<Derived, typename Derived::sub_functor>(pos_tuple, sub_caches, this_copy.argumentsTuple));

		retval.argumentsTuple = this_copy.argumentsTuple;
		retval.trim();

		return retval;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<std::vector<Derived> > NamedSeries<PIRANHA_NAMED_SERIES_TP>::split(const int n) const
	{
		if (n < 0 || n >= boost::tuples::length<ArgsTupleType>::value) 
		{
			PIRANHA_THROW(value_error,"splitting level must be a non-negative integer less than the echelon level of the series");
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


	/// Flatten series.
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
		const std::vector<typename Derived::TermType> tmp(derivedConstCast->flattenTerms(argumentsTuple));
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
