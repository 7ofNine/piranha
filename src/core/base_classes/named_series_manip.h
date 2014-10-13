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
#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

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
	struct named_series_swap {
		
        static void run(ArgsTuple &a1, ArgsTuple &a2)
        {
			a1.get_head().swap(a2.get_head());
			named_series_swap<typename ArgsTuple::tail_type>::run(a1.get_tail(), a2.get_tail());
		}
	};


	template <>
	struct named_series_swap<boost::tuples::null_type> {

		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	/// Swap contents of series.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::swap(Derived &ps2)
	{
		// Do something only if we are not swapping with self.
		if (derivedCast != &ps2)
        {
			named_series_swap<ArgsTupleType>::run(argumentsTuple, ps2.argumentsTuple);
			derivedCast->base_swap(ps2);
		}
	}


	// Meta-programming for appending an argument.
	template <class ArgsDescr>
	struct named_series_append_arg {

		static void run(const std::string &s, typename Ntuple<VectorPsym, boost::tuples::length<ArgsDescr>::value>::type &argsTuple,
						const Psym &arg)
		{
			if (ArgsDescr::head_type::name == s) 
			{
				// Check that the argument is not already present in this set.
				for (VectorPsym::iterator it = argsTuple.get_head().begin(); it != argsTuple.get_head().end(); ++it) 
				{
					if (arg == (*it)) 
					{
						std::cout << "Error: " << s << " argument '" << it->get_name() << "' already present in the set.\n";
						return;
					}
				}
				argsTuple.get_head().push_back(arg);

			} else 
			{
				named_series_append_arg<typename ArgsDescr::tail_type>::run(s, argsTuple.get_tail(), arg);
			}
		}
	};


	template <>
	struct named_series_append_arg<boost::tuples::null_type> {

		static void run(const std::string &s, const boost::tuples::null_type &, const Psym &) 
		{
			std::cout << "Error: '" << s << "' arguments are not known." << std::endl;
		}
	};


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::append_arg(const std::string &s, const Psym &arg)
	{
		PIRANHA_ASSERT(derivedConstCast->empty());
		named_series_append_arg<arguments_description>::run(s, argumentsTuple, arg);
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::append_arg(const Psym &arg)
	{
		PIRANHA_STATIC_CHECK(N >= 0, "Trying to append argument with invalid index.");
		PIRANHA_ASSERT(derivedConstCast->empty());
		// Check that the argument is not already present in this set.
		for (VectorPsym::iterator it = argumentsTuple.template get<N>().begin(); it != argumentsTuple.template get<N>().end(); ++it) 
		{
			if (arg == (*it)) 
			{
				std::cout << "Error: argument '" << it->get_name() << "' already present in the set." << std::endl;
				return;
			}
		}

		argumentsTuple.template get<N>().push_back(arg);
	}



	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::merge_args(const Derived2 &ps2)
	{
		if (static_cast<void *>(this) == static_cast<void const *>(&ps2)) 
		{
			__PDEBUG(std::cout << "Trying to merge with self, returning." << std::endl);
			return;
		}

		if (unlikely(!is_args_compatible(ps2))) 
		{
			merge_incompatible_args(ps2);
		}
	}


	// Template metaprogramming for getting the relative layout of two series.
	// used during multiplication. shomehow generates the knowledge wich parameter is where 
	template <class ArgsTuple>
	struct named_series_get_layout {

		static void run(const ArgsTuple &a1, const ArgsTuple &a2, typename Ntuple < std::vector<std::pair<bool, std::size_t> >,
						boost::tuples::length<ArgsTuple>::value >::type &layout) 
		{
			// Store frequently-used variables.
			const VectorPsym &v1   = a1.get_head(); 
			const VectorPsym &v2   = a2.get_head();
			const std::size_t size1 = v1.size(); 
			const std::size_t size2 = v2.size();
			std::vector<std::pair<bool, std::size_t> > &l = layout.get_head();
			// First we must construct v2's layout wrt to v1.
			l.resize(size2);
			for (std::size_t i = 0; i < size2; ++i) 
            {
				// Let's see if current v2's symbol is present in v1.
				const VectorPsym::const_iterator result = std::find(v1.begin(), v1.end(), v2[i]);
				if (result == v1.end()) 
                {
					l[i].first = false;
				} else 
                {
					// If present, mark its position in v1
					l[i].first = true;
					l[i].second = result - v1.begin();
				}
			}

			// Now we must take care of those elements of v1 that are not represented in
			// the layout (i.e., they are not in v2)
			for (std::size_t i = 0; i < size1; ++i) 
            {
				// Look for element index i in the layout.
				bool found = false;
				const std::size_t l_size = l.size();
				for (std::size_t j = 0; j < l_size; ++j) 
                {
					if (l[j].first && l[j].second == i) 
                    {
						found = true;
						break;
					}
				}

				// If we did not find it, append it to the layout.
				if (!found) 
                {
					l.push_back(std::pair<bool, std::size_t>(true, i));
				}
			}

			named_series_get_layout<typename ArgsTuple::tail_type>::run(a1.get_tail(), a2.get_tail(), layout.get_tail());
		}
	};


	template <>
	class named_series_get_layout<boost::tuples::null_type>
	{
		public:
			static void run(const boost::tuples::null_type &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};



	// Template metaprogramming for applying a layout to a series.
	template <class ArgsTuple>
	struct named_series_apply_layout_to_args {

		static void run(ArgsTuple &a1, const ArgsTuple &a2, const typename Ntuple < std::vector<std::pair<bool, std::size_t> >,
						boost::tuples::length<ArgsTuple>::value >::type &layout) 
		{
			// Store frequently-used variables.
			VectorPsym       &v1 = a1.get_head();
			const VectorPsym &v2 = a2.get_head();
			const std::vector<std::pair<bool, std::size_t> > &l = layout.get_head();
			const std::size_t l_size = l.size();
			// The layout must have at least all arguments in v1.
			PIRANHA_ASSERT(l_size >= v1.size());
			// Memorize the old vector.
			const VectorPsym old(v1);
			// Make space.
			v1.reserve(l_size);

			for (std::size_t i = 0; i < l_size; ++i) 
            {
				if (l[i].first) 
                {
					// The argument was present in the old arguments sets. Copy it over.
					PIRANHA_ASSERT(l[i].second < old.size());
					if (i < v1.size())
					{
						v1[i] = old[l[i].second];

					} else 
					{
						v1.push_back(old[l[i].second]);
					}
				} else 
                {
					// The argument was not present in the old arguments sets. Fetch it from a2.
					PIRANHA_ASSERT(i < v2.size());
					if (i < v1.size()) 
                    {
						v1[i] = v2[i];

					} else 
                    {
						v1.push_back(v2[i]);
					}
				}
			}

			named_series_apply_layout_to_args<typename ArgsTuple::tail_type>::run(a1.get_tail(), a2.get_tail(), layout.get_tail());
		}
	};


	template <>
	class named_series_apply_layout_to_args<boost::tuples::null_type>
	{
		public:
			static void run(const boost::tuples::null_type &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::merge_incompatible_args(const Derived2 &ps2)
	{
		// Build an empty retval and assign it the same arguments as this.
		Derived retval;
		retval.argumentsTuple = argumentsTuple;
		// Build a tuple of layouts.
		typename Ntuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelon_level + 1>::type l;

		// Get the relative layouts of this wrt ps2 and put the result into l.
		named_series_get_layout<ArgsTupleType>::run(retval.argumentsTuple, ps2.arguments(), l);
		
		// Apply the layout to the arguments tuple of retval.
		named_series_apply_layout_to_args<ArgsTupleType>::run(retval.argumentsTuple, ps2.arguments(), l);
		
		// Apply the layout to all terms of this, which will be inserted into retval.
		derivedConstCast->apply_layout_to_terms(l, retval, retval.argumentsTuple);
		
		// Finally, swap the contents of retval with this.
		swap(retval);
	}


	template <class TrimFlags, class ArgsTuple>
	struct trim_flags_init {

		static void run(TrimFlags &tf, const ArgsTuple &argsTuple)
        {
			const std::size_t size = argsTuple.get_head().size();
			tf.get_head().resize(size);
			for (std::size_t i = 0; i < size; ++i) 
            {
				tf.get_head()[i] = false;
			}

			trim_flags_init<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(tf.get_tail(), argsTuple.get_tail());
		}
	};


	template <>
	struct trim_flags_init<boost::tuples::null_type, boost::tuples::null_type> {

		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	inline bool trim_flags_proceed(const boost::tuples::null_type &)
	{
		return false;
	}


	template <class TrimFlags>
	inline bool trim_flags_proceed(const TrimFlags &tf)
	{
		const std::size_t size = tf.get_head().size();
		for (std::size_t i = 0; i < size; ++i) 
		{
			// If we find a flag that was never turned on, we have something to trim.
			if (!tf.get_head()[i]) 
			{
				return true;
			}
		}
		return trim_flags_proceed(tf.get_tail());
	}


	template <class TrimFlags, class ArgsTuple>
	struct trim_arguments 
	{
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
			trim_arguments<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(tf.get_tail(), argsTuple.get_tail());
		}
	};


	template <>
	struct trim_arguments<boost::tuples::null_type, boost::tuples::null_type> 
	{
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::trim()
	{
		typedef typename Ntuple<std::vector<char>, Derived::echelon_level + 1>::type trim_flags_type;
		trim_flags_type trim_flags;
		trim_flags_init<trim_flags_type, ArgsTupleType>::run(trim_flags, argumentsTuple);
		derivedConstCast->trim_test_terms(trim_flags);

		if (trim_flags_proceed(trim_flags)) 
		{
			// First let's do the arguments.
			trim_arguments<trim_flags_type, ArgsTupleType>::run(trim_flags, argumentsTuple);
			// Let's proceed to the terms now.
			Derived tmp;
			derivedCast->trim_terms(trim_flags, tmp, argumentsTuple);
			derivedCast->base_swap(tmp);
		}
	}


	template <class SubCaches, class SubSeries, class ArgsTuple>
	struct init_sub_caches
	{
		static void run(SubCaches &sub_caches, const SubSeries &s, const ArgsTuple *argsTuple) 
		{
			sub_caches.get_head().setup(s, argsTuple);
			init_sub_caches<typename SubCaches::tail_type, SubSeries,ArgsTuple>::run(sub_caches.get_tail(), s, argsTuple);
		}
	};


	template <class SubSeries, class ArgsTuple>
	struct init_sub_caches<boost::tuples::null_type,SubSeries, ArgsTuple>
	{
		static void run(const boost::tuples::null_type &, const SubSeries &, const ArgsTuple *) {}
	};

    //
    // substitute series s for argument name name
    //
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class SubSeries>
	inline Derived NamedSeries<__PIRANHA_NAMED_SERIES_TP>::sub(const std::string &name, const SubSeries &s) const
	{
		typedef typename Derived::term_type::cf_type::
			template sub_cache_selector<SubSeries, typename Derived::term_type::key_type::
			template sub_cache_selector<SubSeries, boost::tuples::null_type, ArgsTupleType>
			::type, ArgsTupleType>::type    sub_caches_type;

		typedef typename Ntuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelon_level + 1>::type    pos_tuple_type;

		PIRANHA_STATIC_CHECK(boost::tuples::length<sub_caches_type>::value == boost::tuples::length<pos_tuple_type>::value,
			"Size mismatch for position and cache tuples in series substitution.");

		const Psym p(name);
		sub_caches_type sub_caches;
		Derived this_copy(*derivedConstCast);
		SubSeries s_copy(s);
		this_copy.merge_args(s_copy);
		s_copy.merge_args(this_copy);

		// Init sub caches using s_copy and this_copy.m_arguments.
		init_sub_caches<sub_caches_type, SubSeries, ArgsTupleType>::run(sub_caches, s_copy, &this_copy.argumentsTuple);

		const pos_tuple_type pos_tuple = psyms2pos(VectorPsym(1,p), this_copy.argumentsTuple);

		Derived retval(this_copy.template base_sub<Derived, typename Derived::sub_functor>(pos_tuple, sub_caches, this_copy.argumentsTuple));

		retval.argumentsTuple = this_copy.argumentsTuple;
		retval.trim();

		return retval;
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<std::vector<Derived> > NamedSeries<__PIRANHA_NAMED_SERIES_TP>::split(const int &n) const
	{
		if (n < 0 || n >= boost::tuples::length<ArgsTupleType>::value) 
		{
			PIRANHA_THROW(value_error,"splitting level must be a non-negative integer less than the echelon level of the series");
		}

		std::vector<std::vector<Derived> > retval;
		derivedConstCast->base_split(retval, n, argumentsTuple);
		const std::size_t size = retval.size();
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
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<Derived> NamedSeries<__PIRANHA_NAMED_SERIES_TP>::flatten() const
	{
		const std::vector<typename Derived::term_type> tmp(derivedConstCast->flatten_terms(argumentsTuple));
		std::vector<Derived> retval;
		const typename std::vector<typename Derived::term_type>::const_iterator itf(tmp.end());
		for  (typename std::vector<typename Derived::term_type>::const_iterator it = tmp.begin(); it != itf; ++it) 
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
