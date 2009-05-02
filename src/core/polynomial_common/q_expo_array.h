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

#ifndef PIRANHA_Q_EXPO_ARRAY_H
#define PIRANHA_Q_EXPO_ARRAY_H

#include <boost/algorithm/string/split.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <memory> // For default allocator.
#include <string>
#include <vector>

#include "../base_classes/q_array.h"
#include "../base_classes/toolbox.h"
#include "../mp.h"
#include "../utils.h"

#define __PIRANHA_Q_EXPO_ARRAY_TP_DECL int Pos, class Allocator
#define __PIRANHA_Q_EXPO_ARRAY_TP Pos,Allocator

namespace piranha
{
	template < __PIRANHA_Q_EXPO_ARRAY_TP_DECL = std::allocator<char> >
	struct q_expo_array {
		typedef toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > type;
	};

	/// Rational exponents array.
	/**
	 * It wraps a piranha::q_array and adds the
	 * capabilities needed for exponent manipulation.
	 */
	template < __PIRANHA_Q_EXPO_ARRAY_TP_DECL >
	class toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> >: public q_array<__PIRANHA_Q_EXPO_ARRAY_TP, toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > >
	{
			typedef q_array<__PIRANHA_Q_EXPO_ARRAY_TP, toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > > ancestor;
			friend class q_array<__PIRANHA_Q_EXPO_ARRAY_TP, toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > >;
		public:
			typedef typename ancestor::value_type value_type;
			typedef typename ancestor::size_type size_type;
			typedef double eval_type;
// 			template <class SubSeries, class SubCachesCons, class ArgsTuple>
// 			struct sub_cache_selector {
// 				typedef boost::tuples::cons<sub_cache<SubSeries,ArgsTuple>,
// 					SubCachesCons> type;
// 			};
// 			template <class SubSeries, class SubCachesCons, class ArgsTuple>
// 			struct ei_sub_cache_selector {
// 				typedef boost::tuples::cons<ei_sub_cache<SubSeries,ArgsTuple>,
// 					SubCachesCons> type;
// 			};
			// Ctors.
			/// Default ctor.
			toolbox(): ancestor() {}
			/// Ctor from string.
			template <class ArgsTuple>
			explicit toolbox(const std::string &s, const ArgsTuple &): ancestor() {
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, this->separator)));
				const size_type w = boost::numeric_cast<size_type>(sd.size());
				this->resize(w);
				for (size_t i = 0; i < w; ++i) {
					(*this)[i] = utils::lexical_converter<value_type>(sd[i]);
				}
			}
			/// Ctor from psym.
			template <class ArgsTuple>
			explicit toolbox(const psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}
			/// Exponentiation to double.
			template <class ArgsTuple>
			toolbox pow(const double &x, const ArgsTuple &args_tuple) const
			{
				return pow_impl(x,args_tuple);
			}
			/// Exponentiation to rational.
			template <class ArgsTuple>
			toolbox pow(const mp_rational &q, const ArgsTuple &args_tuple) const
			{
				return pow_impl(q,args_tuple);
			}
		private:
			template <class T, class ArgsTuple>
			toolbox pow_impl(const T &x, const ArgsTuple &) const
			{
				toolbox retval(*this);
				if (x == -1) {
					retval.invert_sign();
				} else {
					const size_type w = this->size();
					for (size_type i = 0; i < w; ++i) {
						retval[i] *= x;
					}
				}
				return retval;
			}
	};
}

#undef __PIRANHA_Q_EXPO_ARRAY_TP_DECL
#undef __PIRANHA_Q_EXPO_ARRAY_TP

#endif
