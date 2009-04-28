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
#include <string>
#include <utility>
#include <vector>

#include "../config.h" // For (un)likely
#include "../exceptions.h"
#include "../ntuple.h"
#include "../psym.h"

namespace piranha
{
	// Template metaprogramming for swapping of arguments in named series.
	template <class ArgsTuple>
	struct named_series_swap {
		static void run(ArgsTuple &a1, ArgsTuple &a2) {
			a1.get_head().swap(a2.get_head());
			named_series_swap<typename ArgsTuple::tail_type>::run(a1.get_tail(),
					a2.get_tail());
		}
	};

	template <>
	struct named_series_swap<boost::tuples::null_type> {
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

	/// Swap contents of series.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::swap(Derived &ps2)
	{
		// Do something only if we are not swapping with self.
		if (derived_cast != &ps2) {
			named_series_swap<args_tuple_type>::run(m_arguments, ps2.m_arguments);
			derived_cast->swap_terms(ps2);
		}
	}

	// Meta-programming for appending an argument.
	template <class ArgsDescr>
	struct named_series_append_arg {
		static void run(const std::string &s,
			typename ntuple<vector_psym, boost::tuples::length<ArgsDescr>::value>::type &args_tuple,
			const psym &arg) {
			if (ArgsDescr::head_type::name == s) {
				// Check that the argument is not already present in this set.
				for (vector_psym::iterator it = args_tuple.get_head().begin(); it != args_tuple.get_head().end(); ++it) {
					if (arg == (*it)) {
						std::cout << "Error: " << s << " argument '" << it->get_name() << "' already present in the set.\n";
						return;
					}
				}
				args_tuple.get_head().push_back(arg);
			} else {
				named_series_append_arg<typename ArgsDescr::tail_type>::run(s, args_tuple.get_tail(), arg);
			}
		}
	};

	template <>
	struct named_series_append_arg<boost::tuples::null_type> {
		static void run(const std::string &s, const boost::tuples::null_type &, const psym &) {
			std::cout << "Error: '" << s << "' arguments are not known." << std::endl;
		}
	};

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::append_arg(const std::string &s, const psym &arg)
	{
		p_assert(derived_const_cast->empty());
		named_series_append_arg<arguments_description>::run(s, m_arguments, arg);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::append_arg(const psym &arg)
	{
		p_static_check(N >= 0, "Trying to append argument with invalid index.");
		p_assert(derived_const_cast->empty());
		// Check that the argument is not already present in this set.
		for (vector_psym::iterator it = m_arguments.template get<N>().begin();
			it != m_arguments.template get<N>().end(); ++it) {
			if (arg == (*it)) {
				std::cout << "Error: argument '" << it->get_name() << "' already present in the set." << std::endl;
				return;
			}
		}
		m_arguments.template get<N>().push_back(arg);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::merge_args(const Derived2 &ps2)
	{
		if (static_cast<void *>(this) == static_cast<void const *>(&ps2)) {
			__PDEBUG(std::cout << "Trying to merge with self, returning." << std::endl);
			return;
		}
		if (unlikely(!is_args_compatible(ps2))) {
			merge_incompatible_args(ps2);
		}
	}

	// Template metaprogramming for getting the relative layout of two series.
	template <class ArgsTuple>
	struct named_series_get_layout {
		static void run(const ArgsTuple &a1, const ArgsTuple &a2,
			typename ntuple < std::vector<std::pair<bool, size_t> >,
			boost::tuples::length<ArgsTuple>::value >::type &layout) {
			// Store frequently-used variables.
			const vector_psym &v1 = a1.get_head(), &v2 = a2.get_head();
			const size_t size1 = v1.size(), size2 = v2.size();
			std::vector<std::pair<bool, size_t> > &l = layout.get_head();
			// First we must construct v2's layout wrt to v1.
			l.resize(size2);
			for (size_t i = 0; i < size2; ++i) {
				// Let's see if current v2's symbol is present in v1.
				const vector_psym::const_iterator result = std::find(v1.begin(), v1.end(), v2[i]);
				if (result == v1.end()) {
					l[i].first = false;
				} else {
					// If present, mark its position.
					l[i].first = true;
					l[i].second = result - v1.begin();
				}
			}
			// Now we must take care of those elements of v1 that are not represented in
			// the layout (i.e., they are not in v2)
			for (size_t i = 0; i < size1; ++i) {
				// Look for element index i in the layout.
				bool found = false;
				const size_t l_size = l.size();
				for (size_t j = 0; j < l_size; ++j) {
					if (l[j].first && l[j].second == i) {
						found = true;
						break;
					}
				}
				// If we did not find it, append it to the layout.
				if (!found) {
					l.push_back(std::pair<bool, size_t>(true, i));
				}
			}
			named_series_get_layout<typename ArgsTuple::tail_type>::run(a1.get_tail(),
				a2.get_tail(), layout.get_tail());
		}
	};

	template <>
	class named_series_get_layout<boost::tuples::null_type>
	{
		public:
			static void run(const boost::tuples::null_type &, const boost::tuples::null_type &,
				const boost::tuples::null_type &) {}
	};

	// Template metaprogramming for applying a layout to a series.
	template <class ArgsTuple>
	struct named_series_apply_layout_to_args {
		static void run(ArgsTuple &a1, const ArgsTuple &a2, const typename ntuple < std::vector<std::pair<bool, size_t> >,
			boost::tuples::length<ArgsTuple>::value >::type &layout) {
			// Store frequently-used variables.
			vector_psym &v1 = a1.get_head();
			const vector_psym &v2 = a2.get_head();
			const std::vector<std::pair<bool, size_t> > &l = layout.get_head();
			const size_t l_size = l.size();
			// The layout must have at least all arguments in v1.
			p_assert(l_size >= v1.size());
			// Memorize the old vector.
			const vector_psym old(v1);
			// Make space.
			v1.reserve(l_size);
			for (size_t i = 0; i < l_size; ++i) {
				if (l[i].first) {
					// The argument was present in the old arguments sets. Copy it over.
					p_assert(l[i].second < old.size());
					if (i < v1.size()) {
						v1[i] = old[l[i].second];
					} else {
						v1.push_back(old[l[i].second]);
					}
				} else {
					// The argument was not present in the old arguments sets. Fetch it from a2.
					p_assert(i < v2.size());
					if (i < v1.size()) {
						v1[i] = v2[i];
					} else {
						v1.push_back(v2[i]);
					}
				}
			}
			named_series_apply_layout_to_args<typename ArgsTuple::tail_type>::run(a1.get_tail(),
				a2.get_tail(), layout.get_tail());
		}
	};

	template <>
	class named_series_apply_layout_to_args<boost::tuples::null_type>
	{
		public:
			static void run(const boost::tuples::null_type &, const boost::tuples::null_type &,
				const boost::tuples::null_type &) {}
	};

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline void toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::merge_incompatible_args(const Derived2 &ps2)
	{
		// Build an empty retval and assign it the same arguments as this.
		Derived retval;
		retval.m_arguments = m_arguments;
		// Build a tuple of layouts.
		typename ntuple<std::vector<std::pair<bool, size_t> >, n_arguments_sets>::type l;
		// Get the relative layouts of this wrt ps2 and put the result into l.
		named_series_get_layout<args_tuple_type>::run(retval.m_arguments, ps2.arguments(), l);
		// Apply the layout to the arguments tuple of retval.
		named_series_apply_layout_to_args<args_tuple_type>::run(retval.m_arguments, ps2.arguments(), l);
		// Apply the layout to all terms of this, which will be inserted into retval.
		derived_const_cast->apply_layout_to_terms(retval.m_arguments, l, retval);
		// Finally, swap the contents of retval with this.
		swap(retval);
	}

	template <class TrimFlags, class ArgsTuple>
	struct trim_flags_init {
		static void run(TrimFlags &tf, const ArgsTuple &args_tuple) {
			const size_t size = args_tuple.get_head().size();
			tf.get_head().resize(size);
			for (size_t i = 0; i < size; ++i) {
				tf.get_head()[i] = false;
			}
			trim_flags_init<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(
				tf.get_tail(),
				args_tuple.get_tail()
			);
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
		const size_t size = tf.get_head().size();
		for (size_t i = 0; i < size; ++i) {
			// If we find a flag that was never turned on, we have something to trim.
			if (!tf.get_head()[i]) {
				return true;
			}
		}
		return trim_flags_proceed(tf.get_tail());
	}

	template <class TrimFlags, class ArgsTuple>
	struct trim_arguments {
		static void run(const TrimFlags &tf, ArgsTuple &args_tuple) {
			const size_t size = tf.get_head().size();
			p_assert(size == args_tuple.get_head().size());
			vector_psym new_vector;
			for (size_t i = 0; i < size; ++i) {
				if (tf.get_head()[i]) {
					new_vector.push_back(args_tuple.get_head()[i]);
				}
			}
			new_vector.swap(args_tuple.get_head());
			trim_arguments<typename TrimFlags::tail_type, typename ArgsTuple::tail_type>::run(
				tf.get_tail(), args_tuple.get_tail());
		}
	};

	template <>
	struct trim_arguments<boost::tuples::null_type, boost::tuples::null_type> {
		static void run(const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::trim()
	{
		typedef typename ntuple<std::vector<bool>, n_arguments_sets>::type trim_flags_type;
		trim_flags_type trim_flags;
		trim_flags_init<trim_flags_type, args_tuple_type>::run(trim_flags, m_arguments);
		derived_const_cast->trim_test_terms(trim_flags);
		if (trim_flags_proceed(trim_flags)) {
			// First let's do the arguments.
			trim_arguments<trim_flags_type, args_tuple_type>::run(trim_flags, m_arguments);
			// Let's proceed to the terms now.
			Derived tmp;
			derived_cast->trim_terms(trim_flags, tmp, m_arguments);
			derived_cast->swap_terms(tmp);
		}
	}

// 	template <__PIRANHA_NAMED_SERIES_TP_DECL>
// 	template <class Filter>
// 	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::filter(const Filter &f) const
// 	{
// 		typedef typename Derived::const_iterator const_iterator;
// 		Derived retval;
// 		retval.m_arguments = m_arguments;
// 		shared_args::set(m_arguments);
// 		const const_iterator it_f = derived_const_cast->end();
// 		for (const_iterator it =  derived_const_cast->begin(); it != it_f; ++it) {
// 			if (f(*it)) {
// 				retval.insert(*it, retval.m_arguments);
// 			}
// 		}
// 		retval.trim();
// 		return retval;
// 	}

	template <class SubCaches, class SubSeries, class ArgsTuple>
	struct init_sub_caches
	{
		static void run(SubCaches &sub_caches, const SubSeries &s, const ArgsTuple *args_tuple) {
			sub_caches.get_head().setup(s,args_tuple);
			init_sub_caches<typename SubCaches::tail_type,SubSeries,ArgsTuple>::
				run(sub_caches.get_tail(),s,args_tuple);
		}
	};

	template <class SubSeries, class ArgsTuple>
	struct init_sub_caches<boost::tuples::null_type,SubSeries, ArgsTuple>
	{
		static void run(const boost::tuples::null_type &, const SubSeries &, const ArgsTuple *) {}
	};

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class SubSeries>
	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::sub(const psym &p, const SubSeries &s) const
	{
		typedef typename Derived::term_type::cf_type::
			template sub_cache_selector<SubSeries,typename Derived::term_type::key_type::
			template sub_cache_selector<SubSeries,boost::tuples::null_type,args_tuple_type>
			::type,args_tuple_type>::type sub_caches_type;
		typedef typename ntuple<std::vector<std::pair<bool, size_t> >, n_arguments_sets>::type pos_tuple_type;
		p_static_check(boost::tuples::length<sub_caches_type>::value == boost::tuples::length<pos_tuple_type>::value,
			"Size mismatch for position and cache tuples in series substitution.");
		sub_caches_type sub_caches;
		Derived this_copy(*derived_const_cast);
		SubSeries s_copy(s);
		this_copy.merge_args(s_copy);
		s_copy.merge_args(this_copy);
		// Init sub caches using s_copy and this_copy.m_arguments.
		init_sub_caches<sub_caches_type,SubSeries,args_tuple_type>::run(sub_caches,s_copy,&this_copy.m_arguments);
		const pos_tuple_type pos_tuple = psyms2pos(vector_psym(1,p),this_copy.m_arguments);
		Derived retval(this_copy.template base_sub<Derived,typename Derived::sub_functor>(
			pos_tuple, sub_caches, this_copy.m_arguments));
		retval.m_arguments = this_copy.m_arguments;
		retval.trim();
		return retval;
	}
}

#endif
