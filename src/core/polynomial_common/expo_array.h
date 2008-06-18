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

#ifndef PIRANHA_EXPO_ARRAY_H
#define PIRANHA_EXPO_ARRAY_H

#include <memory> // For standard allocator.
#include <string>
#include <vector>

#include "../base_classes/int_array.h"
#include "../psym.h"
#include "expo_array_commons.h"

#define __PIRANHA_EXPO_ARRAY_TP_DECL int Bits, int Pos, class Allocator
#define __PIRANHA_EXPO_ARRAY_TP Bits,Pos,Allocator

namespace piranha
{
	/// Exponents array, dynamically sized version.
	/**
	 * It wraps a piranha::int_array with signed integer sized Bits, and adds the
	 * capabilities needed for exponent manipulation.
	 */
	template < __PIRANHA_EXPO_ARRAY_TP_DECL = std::allocator<char> >
	class expo_array:
				public int_array<Bits, Pos, Allocator, expo_array<__PIRANHA_EXPO_ARRAY_TP> >,
				public expo_array_commons<expo_array<__PIRANHA_EXPO_ARRAY_TP> >
	{
			friend class expo_array_commons<expo_array<__PIRANHA_EXPO_ARRAY_TP> >;
			typedef expo_array_commons<expo_array<__PIRANHA_EXPO_ARRAY_TP> > expo_commons;
			typedef int_array<Bits, Pos, Allocator, expo_array<__PIRANHA_EXPO_ARRAY_TP> > ancestor;
		public:
			typedef typename ancestor::value_type value_type;
			typedef typename ancestor::size_type size_type;
			class proxy: public ancestor::reference_proxy
			{
					typedef typename ancestor::reference_proxy proxy_ancestor;
				public:
					typedef proxy type;
					proxy(const expo_array &e): proxy_ancestor(e) {}
					// Expo-array specifics.
					void multiply(proxy e2, expo_array &ret) const {
						proxy_ancestor::m_ptr->multiply(e2, ret);
					}
			};
			// Ctors.
			/// Default ctor.
			expo_array(): ancestor::int_array() {}
			/// Ctor from string.
			template <class ArgsTuple>
			explicit expo_array(const std::string &s, const ArgsTuple &): ancestor::int_array(),
					expo_commons::expo_array_commons(s) {}
			/// Ctor from psym.
			template <class ArgsTuple>
			explicit expo_array(const psym_p &p, const int &n, const ArgsTuple &a): ancestor::int_array(p, n, a) {}
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
			template <class ExpoArray, class ResultType>
			void multiply(const ExpoArray &e2, ResultType &ret) const {
				const size_type max_w = ancestor::size(), min_w = e2.size();
				// Resize, if needed.
				ret.resize(max_w);
				// Assert widths, *this should always come from a polynomial, and its width should hence be
				// already adjusted my merge_args in multiplication routines.
				p_assert(max_w >= min_w);
				p_assert(ret.size() == max_w);
				size_type i;
				for (i = 0;i < min_w;++i) {
					ret[i] = (*this)[i] + e2[i];
				}
				for (;i < max_w;++i) {
					ret[i] = (*this)[i];
				}
			}
	};
}

#undef __PIRANHA_EXPO_ARRAY_TP_DECL
#undef __PIRANHA_EXPO_ARRAY_TP

#endif
