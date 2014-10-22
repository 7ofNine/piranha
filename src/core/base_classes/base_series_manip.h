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
	 */
	template <class Term2, class Term1>
	class term_converter
	{
		public:
			/// Constructor.
			template <class ArgsTuple>
			explicit term_converter(const Term2 &c, const ArgsTuple &a): result(c,a) {}
			/// Copy of the converted term.
			const Term1 result;
	};

	/// Specialized term converter.
	/**
	 * It will be invoked when the type to convert from is the same as the converted type. A reference
	 * to the convertee is stored inside the class.
	 */
	template <class Term2>
	class term_converter<Term2, Term2>
	{
		public:
			/// Constructor.
			template <class ArgsTuple>
			explicit term_converter(const Term2 &c, const ArgsTuple &): result(c) {}
			/// Reference to the converted term.
			const Term2 &result;
	};

	// TODO: update doc here.
	/// High-level insertion function.
	/**
	 * This function is used to insert terms into a series. It requires that the number of arguments
	 * of each element of the term is smaller or equal to the series',
	 * otherwise an assertion fails and the program aborts. base_pseries::merge_args,
	 * base_pseries::append_cf_args, base_pseries::append_trig_args, etc. can be used to add the needed arguments
	 * to the series.
	 *
	 * This function performs some checks and then calls llInsert.
	 */
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool CanonicalCheck, bool Sign, class Term2, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::insert(const Term2 &term_, const ArgsTuple &argsTuple)
	{
		term_converter<Term2, TermType> converted_term(term_, argsTuple);
		// Make sure the appropriate routines for the management of arguments have been called.
		PIRANHA_ASSERT(converted_term.result.cf.is_insertable(argsTuple) && converted_term.result.key.is_insertable(argsTuple));

		TermType *new_term(0);
		if (unlikely(converted_term.result.cf.needs_padding(argsTuple) ||
			converted_term.result.key.needs_padding(argsTuple))) 
		{
			new_term = TermType::allocator.allocate(1);
			TermType::allocator.construct(new_term, converted_term.result);
			new_term->cf.pad_right(argsTuple);
			new_term->key.pad_right(argsTuple);
		}

		if (CanonicalCheck) 
		{
			if (!converted_term.result.is_canonical(argsTuple)) 
			{
				if (new_term == 0) 
				{
					new_term = TermType::allocator.allocate(1);
					TermType::allocator.construct(new_term, converted_term.result);
				}
				new_term->canonicalise(argsTuple);
			}
		}
		const TermType *insert_term(0);
		if (new_term) 
		{
			insert_term = new_term;
		} else 
		{
			insert_term = &converted_term.result;
		}
		llInsert<Sign>(*insert_term, argsTuple);
		
		if (new_term) 
		{
			TermType::allocator.destroy(new_term);
			TermType::allocator.deallocate(new_term, 1);
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
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::insertRange(const Iterator &begin,
		const Iterator &end, const ArgsTuple &argsTuple)
	{
		for (Iterator it = begin; it != end; ++it) 
		{
			insert(*it,argsTuple);
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	inline typename BaseSeries<__PIRANHA_BASE_SERIES_TP>::const_iterator
	BaseSeries<__PIRANHA_BASE_SERIES_TP>::findTerm(const TermType &t) const
	{
		return container.find(t);
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <bool Sign, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::llInsert(const TermType &term, const ArgsTuple &argsTuple)
	{
		// TODO: think about moving this check higher in the stack of functions, we probably don't want to reach
		// _this_ point before checking for ignorability.
		if (term.cf.is_ignorable(argsTuple) || term.key.is_ignorable(argsTuple)) 
		{
			return;
		}
		PIRANHA_ASSERT(term.cf.is_insertable(argsTuple) && term.key.is_insertable(argsTuple) &&
			!term.cf.needs_padding(argsTuple) && !term.key.needs_padding(argsTuple) && term.is_canonical(argsTuple));

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
			if (it->cf.is_ignorable(argsTuple)) 
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
		std::pair<const_iterator, bool> res(container.insert(term));
		PIRANHA_ASSERT(res.second);
		if (!Sign) 
		{
			res.first->cf.invert_sign(argsTuple);
		}
	}


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
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSwap(Derived &ps2)
	{
		PIRANHA_ASSERT(derived_cast != &ps2);
		container.swap(ps2.container);
	}


	/// Apply an arguments layout to all terms and insert them into retval.
	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class Layout, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::applyLayoutToTerms(const Layout &l, Derived &retval, const ArgsTuple &argsTuple) const
	{
		const const_iterator itf = end();
		for (const_iterator it = begin(); it != itf; ++it) 
		{
			TermType term(*it);
			term.cf.apply_layout(l, argsTuple);
			term.key.apply_layout(l, argsTuple);
			retval.insert(term, argsTuple);
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::trimTestTerms(TrimFlags &tf) const
	{
		const const_iterator itf = end();
		for (const_iterator it = begin(); it != itf; ++it) 
		{
			it->cf.trim_test(tf);
			it->key.trim_test(tf);
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class TrimFlags, class ArgsTuple>
	inline void BaseSeries<__PIRANHA_BASE_SERIES_TP>::trimTerms(const TrimFlags &tf, Derived &retval, const ArgsTuple &argsTuple) const
	{
		const const_iterator itf = end();
		for (const_iterator it = begin(); it != itf; ++it) 
		{
			retval.insert(TermType(typename TermType::CfType(it->cf.trim(tf, argsTuple)), typename TermType::KeyType(it->key.trim(tf, argsTuple))), argsTuple);
		}
	}


	template <__PIRANHA_BASE_SERIES_TP_DECL>
	template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
	inline RetSeries BaseSeries<__PIRANHA_BASE_SERIES_TP>::baseSub(const PosTuple &pos_tuple, SubCaches &sub_caches, const ArgsTuple &argsTuple) const
	{
		PIRANHA_STATIC_CHECK((boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value), "Positional and arguments' tuples' lengths do not match.");

		RetSeries retval;
		const const_iterator itf = end();
		for (const_iterator it = begin(); it != itf; ++it) 
		{
			RetSeries tmp = SubFunctor::template run<RetSeries>(it->cf, pos_tuple, sub_caches,argsTuple);
			// NOTICE: series multadd here?
			tmp.baseMultBy(SubFunctor::template run<RetSeries>(it->key, pos_tuple, sub_caches, argsTuple), argsTuple);

			retval.baseAdd(tmp,argsTuple);
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
				PIRANHA_THROW(value_error,"cannot split up to the specified level: series is non-degenerate");
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
