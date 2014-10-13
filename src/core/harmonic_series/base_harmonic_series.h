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

#ifndef PIRANHA_BASE_HARMONIC_SERIES_H
#define PIRANHA_BASE_HARMONIC_SERIES_H

#include <algorithm> // For max_element and min_element.

#include "../config.h"
#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

// TODO: share implementation with power series class, using templates and functors?

namespace piranha
{
	/// Base harmonic series toolbox.
	template <int HarmonicArgsPosition, int HarmonicTermPosition, class HarmonicDegree, class Derived>
	class BaseHarmonicSeries
	{
			PIRANHA_STATIC_CHECK(HarmonicArgsPosition >= 0, "Invalid harmonic args position.");

			template <class Term>
			struct HarmonicDegreeBinaryPredicate
			{
				bool operator()(const Term &t1, const Term &t2) const 
				{
					return (t1.template get<HarmonicTermPosition>().harmonicDegree() < t2.template get<HarmonicTermPosition>().harmonicDegree());
				}
			};


			template <class Term>
			struct HarmonicOrderBinaryPredicate
			{
				bool operator()(const Term &t1, const Term &t2) const 
				{
					return (t1.template get<HarmonicTermPosition>().harmonicOrder() < t2.template get<HarmonicTermPosition>().harmonicOrder());
				}
			};


			template <class Term, class PosTuple>
			struct PartialHarmonicDegreeBinaryPredicate
			{
				PartialHarmonicDegreeBinaryPredicate(const PosTuple &p):posTuple(p) {}


				bool operator()(const Term &t1, const Term &t2) const 
				{
					return (t1.template get<HarmonicTermPosition>().partial_h_degree(posTuple) < t2.template get<HarmonicTermPosition>().partial_h_degree(posTuple));
				}

				const PosTuple &posTuple;
			};


			template <class Term, class PosTuple>
			struct partial_h_order_binary_predicate
			{
				partial_h_order_binary_predicate(const PosTuple &p):m_p(p) {}

				bool operator()(const Term &t1, const Term &t2) const 
				{
					return (t1.template get<HarmonicTermPosition>().partial_h_order(m_p) < t2.template get<HarmonicTermPosition>().partial_h_order(m_p));
				}

				const PosTuple &m_p;
			};


		public:

			static const int harmonic_args_position = HarmonicArgsPosition;
			static const int harmonic_term_position = HarmonicTermPosition;
			typedef HarmonicDegree h_degree_type;
			
			/// Get the harmonic degree.
			HarmonicDegree harmonicDegree() const 
			{
				if (derived_const_cast->empty()) 
				{
					return HarmonicDegree(0);
				}
				const typename Derived::const_iterator result(std::max_element(derived_const_cast->begin(), derived_const_cast->end(),
																               HarmonicDegreeBinaryPredicate<typename Derived::term_type>()));

				return result->template get<HarmonicTermPosition>().harmonicDegree();
			}


			/// Get the harmonic order.
			/**
			 * The harmonic order is defined as the minimum harmonic degree of the terms composing the series.
			 */
			HarmonicDegree harmonicOrder() const 
			{
				if (derived_const_cast->empty()) 
				{
					return HarmonicDegree(0);
				}
				const typename Derived::const_iterator result(std::min_element(
							derived_const_cast->begin(),
							derived_const_cast->end(),
							HarmonicOrderBinaryPredicate<typename Derived::term_type>()
						));

				return result->template get<HarmonicTermPosition>().harmonicOrder();
			}


			/// Return true if all terms of the harmonic series are cosines.
			/**
			 * An empty series will return false.
			 */
			bool is_cosine() const
			{
				typedef typename Derived::const_iterator const_iterator;

				if (derived_const_cast->empty()) 
				{
					return false;
				}
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
				{
					if (!it->template get<HarmonicTermPosition>().getFlavour()) 
					{
						return false;
					}
				}

				return true;
			}


			/// Return true if all terms of the harmonic series are sines.
			/**
			 * An empty series will return false.
			 */
			bool is_sine() const
			{
				typedef typename Derived::const_iterator const_iterator;
				if (derived_const_cast->empty()) 
				{
					return false;
				}
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
				{
					if (it->template get<HarmonicTermPosition>().getFlavour()) 
					{
						return false;
					}
				}
				return true;
			}


		//protected:
			/// Get the harmonic degree of the series for specific variables.
			template <class PosTuple>
			HarmonicDegree base_partial_h_degree(const PosTuple &pos_tuple) const 
			{
				if (derived_const_cast->empty()) 
				{
					return HarmonicDegree(0);
				}
				const typename Derived::const_iterator result(std::max_element(
							derived_const_cast->begin(),
							derived_const_cast->end(),
							PartialHarmonicDegreeBinaryPredicate<typename Derived::term_type,PosTuple>(pos_tuple)
						));
				
				return result->template get<HarmonicTermPosition>().partial_h_degree(pos_tuple);
			}


			/// Get the mininum harmonic degree of the series for specific variables.
			template <class PosTuple>
			HarmonicDegree base_partial_h_order(const PosTuple &pos_tuple) const 
			{
				if (derived_const_cast->empty()) 
				{
					return HarmonicDegree(0);
				}

				const typename Derived::const_iterator result(std::min_element(
							                                   derived_const_cast->begin(),
							                                   derived_const_cast->end(),
							                                   partial_h_order_binary_predicate<typename Derived::term_type, PosTuple>(pos_tuple) ));
				
				return result->template get<HarmonicTermPosition>().partial_h_order(pos_tuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
