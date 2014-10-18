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

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <cstddef>
#include <utility>
#include <vector>

#include "../config.h" // For (un)likely
#include "../mp.h"
#include "../Psym.h"
#include "named_series_def.h"
#include "named_series_mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	// TMP for arguments compatibility check.
	template <class ArgsTuple>
	inline bool named_series_is_args_compatible(const ArgsTuple & a1, const ArgsTuple & a2)
	{
		const std::size_t w = a2.get_head().size();

		if (unlikely(w > a1.get_head().size())) 
		{
			return false;
		}

		for (std::size_t i = 0; i < w; ++i) 
		{
			if (unlikely(a1.get_head()[i] != a2.get_head()[i])) 
			{
				return false;
			}
		}

		return named_series_is_args_compatible(a1.get_tail(), a2.get_tail());
	}


	template <>
	inline bool named_series_is_args_compatible<boost::tuples::null_type>(const boost::tuples::null_type &, const boost::tuples::null_type &)
	{
		return true;
	}


	/// Compatibility check for arguments.
	/**
	 * Test whether series' arguments are compatible with those from ps2. Compatibility
	 * means that the number of arguments in all arguments sets are equal to or greater than ps2's, and
	 * that arguments have the same positions as in ps2's.
	 * @param[in] ps2 series compatibility is tested against.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline bool NamedSeries<__PIRANHA_NAMED_SERIES_TP>::is_args_compatible(const Derived2 &ps2) const
	{
		// Use getter in second place because we may be interacting with other series type.
		return named_series_is_args_compatible(argumentsTuple, ps2.arguments());
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline double NamedSeries<__PIRANHA_NAMED_SERIES_TP>::norm() const
	{
		return derived_const_cast->base_norm(argumentsTuple);
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename TermEvalTypeDeterminer<Term>::type
	NamedSeries<__PIRANHA_NAMED_SERIES_TP>::eval(const double &t) const
	{
		return derived_const_cast->base_eval(t, argumentsTuple);
	}

	// TMP function for checking that evaluation dictionary has all the elements needed.
	static inline bool check_eval_dict(const eval_dict &, const boost::tuples::null_type &)
	{
		return true;
	}

	template <class ArgsTuple>
	static inline bool check_eval_dict(const eval_dict &d, const ArgsTuple &argsTuple)
	{
		const std::size_t size = argsTuple.get_head().size();
		const eval_dict::const_iterator it_f = d.end();

		for (std::size_t i = 0; i < size; ++i)
        {
			// If the dictionary does not contain the symbol's name, return false.
			if (d.find(argsTuple.get_head()[i].get_name()) == it_f)
            {
				return false;
			}
		}

		// Check next tuple position.
		return check_eval_dict(d, argsTuple.get_tail());
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename TermEvalTypeDeterminer<Term>::type
	NamedSeries<__PIRANHA_NAMED_SERIES_TP>::eval(const eval_dict &d) const
	{
		if (!check_eval_dict(d, argumentsTuple)) 
		{
			PIRANHA_THROW(value_error,"evaluation dictionary does not contain entries for all the symbols of the series");
		}

		// Vector of original time evaluation vectors, paired with the corresponding symbol.
		// NOTE: here we allocate dynamically. This can be avoided by fixing a max number of items in time evaluation for psyms
		// and using static vectors. We should test performance before bothering though.
		std::vector<std::pair<Psym,std::vector<double> > > orig_eval;
		orig_eval.reserve(d.size());
		const eval_dict::const_iterator it_f(d.end());

		for (eval_dict::const_iterator it = d.begin(); it != it_f; ++it)
        {
			orig_eval.push_back(std::make_pair(Psym(it->first), Psym(it->first).get_time_eval()));
			Psym(it->first, it->second);
		}

		const typename TermEvalTypeDeterminer<Term>::type retval(eval(0));

		// Restore original evaluation vectors.
		// NOTE: here RAII here, to be exception safe?
		for (std::size_t i = 0; i < orig_eval.size(); ++i)
        {
			orig_eval[i].first.set_time_eval(orig_eval[i].second);
		}
		return retval;
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::size_t NamedSeries<__PIRANHA_NAMED_SERIES_TP>::psi(const int &start, const int &step) const
	{
		return derived_const_cast->psi_(start, step, argumentsTuple);
	}


	inline bool tuple_vector_same_sizes(const boost::tuples::null_type &, const boost::tuples::null_type &)
	{
		return true;
	}


	template <class TupleVector1, class TupleVector2>
	inline bool tuple_vector_same_sizes(const TupleVector1 &tv1, const TupleVector2 &tv2)
	{
		if (tv1.get_head().size() != tv2.get_head().size()) 
        {
			return false;
		}
		return tuple_vector_same_sizes(tv1.get_tail(), tv2.get_tail());
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline bool NamedSeries<__PIRANHA_NAMED_SERIES_TP>::series_comparison(const T &other) const
	{
		// If the sizes of the arguments tuples do not coincide, series are different.
		if (!tuple_vector_same_sizes(argumentsTuple, other.argumentsTuple)) 
        {
			return false;
		}

		// If arguments tuples are completely identical, just run the base
		// comparison function.
		if (argumentsTuple == other.arguments()) 
        {
			return derived_const_cast->base_equal_to(other);
		}

		// If we have same sizes of arguments tuples but they are not identical, then we may have to do
		// an arguments merge an see if the arguments are permutated or they are really different.
		// NOTE: this check is repeated in base_equal_to, but doing it here could save a lot of work below.
		if (derived_const_cast->length() != other.length())
        {
			return false;
		}


		// Build a tuple of layouts.
		typename Ntuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::type l;
		// Get the relative layouts of this wrt other and put the result into l.
		named_series_get_layout<ArgsTupleType>::run(argumentsTuple, other.arguments(), l);

		// If the layout is bigger than the current ags tuple, it means that it is not a permutation,
		// there are different arguments in this and other. Hence we can return false.
		if (!tuple_vector_same_sizes(argumentsTuple, l)) 
        {
			return false;
		}

		// In this last case, the arguments are the same but they are ordered differently. Need to copy
		// this into new series with correct ordering and then do the comparison.
		// Build an empty retval and assign it the same arguments as this.
		Derived tmp;
		tmp.argumentsTuple = argumentsTuple;

		// Apply the layout to the arguments tuple of retval.
		named_series_apply_layout_to_args<ArgsTupleType>::run(tmp.argumentsTuple, other.arguments(), l);
		
        // Apply the layout to all terms of this and insert them into tmp.
		derived_const_cast->applyLayoutToTerms(l, tmp, tmp.argumentsTuple);
		
        // Now we can perform the comparison between tmp and other.
		return tmp.base_equal_to(other);
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline bool NamedSeries<__PIRANHA_NAMED_SERIES_TP>::operator==(const T &x) const
	{
		return NamedSeriesEqualitySelector<T>::run(*derived_const_cast, x);
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class T>
	inline bool NamedSeries<__PIRANHA_NAMED_SERIES_TP>::operator!=(const T &x) const
	{
		return !(operator==(x));
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
