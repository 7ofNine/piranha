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

#ifndef PIRANHA_TRIG_ARRAY_H
#define PIRANHA_TRIG_ARRAY_H

#include <boost/tuple/tuple.hpp> // For sub cache selection.
#include <complex> // For std::complex<SubSeries>.
#include <memory> // For standard allocator.
#include <string>

#include "../base_classes/int_array.h"
#include "../common_functors.h"
#include "../int_power_cache.h"
#include "trig_array_commons.h"

#define __PIRANHA_TRIG_ARRAY_TP_DECL int Bits, int Pos, class Allocator
#define __PIRANHA_TRIG_ARRAY_TP Bits,Pos,Allocator

namespace piranha
{
	/// Trigonometric array, dynamically sized version.
	/**
	 * It wraps a piranha::int_array with signed integer sized Bits, and adds the
	 * capabilities needed for trigonometric manipulation.
	 */
	template < __PIRANHA_TRIG_ARRAY_TP_DECL = std::allocator<char> >
	class trig_array:
				public int_array<Bits, Pos, Allocator, trig_array<__PIRANHA_TRIG_ARRAY_TP> >,
				public trig_array_commons<trig_array<__PIRANHA_TRIG_ARRAY_TP> >
	{
			friend class trig_array_commons<trig_array<__PIRANHA_TRIG_ARRAY_TP> >;
			typedef trig_array_commons<trig_array<__PIRANHA_TRIG_ARRAY_TP> > trig_commons;
			typedef int_array<Bits, Pos, Allocator, trig_array<__PIRANHA_TRIG_ARRAY_TP> > ancestor;
		public:
			typedef typename ancestor::value_type value_type;
			typedef typename ancestor::size_type size_type;
			typedef double eval_type;
			class proxy: public ancestor::reference_proxy
			{
					typedef typename ancestor::reference_proxy proxy_ancestor;
					friend class trig_array;
				public:
					proxy(const trig_array &t): proxy_ancestor(t) {}
					// These are trig_array specific.
					bool flavour() const {
						return proxy_ancestor::m_ptr->flavour();
					}
					void multiply(const proxy &t2, trig_array &ret1, trig_array &ret2) const {
						proxy_ancestor::m_ptr->multiply(t2, ret1, ret2);
					}
					bool operator<(const proxy &t2) const {
						return proxy_ancestor::m_ptr->operator<(*t2.m_ptr);
					}
			};
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct sub_cache_selector {
				typedef boost::tuples::cons<int_power_cache<std::complex<SubSeries>,
					base_series_arithmetics<std::complex<SubSeries>,ArgsTuple> >,SubCachesCons> type;
			};
			// Ctors.
			/// Default ctor.
			trig_array(): ancestor::int_array() {}
			/// Ctor from string.
			template <class ArgsTuple>
			explicit trig_array(const std::string &s, const ArgsTuple &): ancestor::int_array(),
					trig_commons::trig_array_commons(s) {}
			template <class ArgsTuple>
			explicit trig_array(const psym_p &p, const int &n, const ArgsTuple &a): ancestor::int_array(p, n, a) {}
			explicit trig_array(const proxy &p) {
				*this = *p.m_ptr;
			}
			template <int Pos2>
			explicit trig_array(const trig_array<Bits,Pos2,Allocator> &ta): ancestor::int_array(ta) {}
			// Probing.
			/// Data footprint.
			/**
			 * Returns the memory occupied by the data members.
			 */
			size_t data_footprint() const {
				return (ancestor::size()*sizeof(value_type));
			}
			// Math.
			/// Multiplication.
			/**
			 * Used in poisson_series_term multiplication.
			 * TODO: update docs below.
			 * Multiplication of two trigonometric functions using Werner's formulas, i.e.
			 * \f[
			 * C\cos\alpha\cdot\cos\beta=
			 * \frac{C}{2} \cos \left( \alpha - \beta \right) + \frac{C}{2} \cos \left( \alpha + \beta \right)
			 * \f]
			 * and the likes. Notice that in the first return value always goes the \f$ \alpha - \beta \f$ term
			 * and in the second one always goes \f$ \alpha + \beta \f$ one.
			 * Please also note that no assumptions are made with respect to return values' content
			 * (e.g., it is not guaranteed that return values are empty).
			 * @param[in] t2 factor.
			 * @param[out] ret1 first return value.
			 * @param[out] ret2 second return value.
			 */
			template <class TrigArray>
			void multiply(const TrigArray &t2, trig_array &ret1, trig_array &ret2) const
			// NOTE: we are not using here a general version of vector addition/subtraction
			// because this way we can do two operations (+ and -) every cycle. This is a performance
			// critical part, so the optimization should be worth the hassle.
			{
				const size_type max_w = ancestor::size(), min_w = t2.size();
				// Assert widths, *this should always come from a regular Poisson series, and its width should hence be
				// already adjusted my merge_args in multiplication routines.
				p_assert(max_w >= min_w);
				// Adjust the width of retvals, if needed.
				ret1.resize(max_w);
				ret2.resize(max_w);
				p_assert(ret1.size() == max_w);
				p_assert(ret2.size() == max_w);
				size_type i;
				for (i = 0; i < min_w; ++i) {
					ret1[i] = (*this)[i] - t2[i];
					ret2[i] = (*this)[i] + t2[i];
				}
				for (; i < max_w; ++i) {
					ret1[i] = (*this)[i];
					ret2[i] = (*this)[i];
				}
			}
	};
}

#undef __PIRANHA_TRIG_ARRAY_TP_DECL
#undef __PIRANHA_TRIG_ARRAY_TP

#endif
