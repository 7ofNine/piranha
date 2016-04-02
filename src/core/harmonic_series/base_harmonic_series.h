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

#define DerivedConstCast static_cast<Derived const *>(this)
//#define derived_cast       static_cast<Derived *>(this)

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
					return (t1.template get<HarmonicTermPosition>().partialHarmonicDegree(posTuple) < t2.template get<HarmonicTermPosition>().partialHarmonicDegree(posTuple));
				}

				const PosTuple &posTuple;
			};


			template <class Term, class PosTuple>
			struct PartialHarmonicOrderBinaryPredicate
			{
				PartialHarmonicOrderBinaryPredicate(const PosTuple &p):posTuple(p) {}

				bool operator()(const Term &t1, const Term &t2) const 
				{
					return (t1.template get<HarmonicTermPosition>().partialHarmonicOrder(posTuple) < t2.template get<HarmonicTermPosition>().partialHarmonicOrder(posTuple));
				}

				const PosTuple &posTuple;
			};


		public:

			static const int harmonicArgsPosition = HarmonicArgsPosition;
			static const int harmonicTermPosition = HarmonicTermPosition;
			typedef HarmonicDegree HarmonicDegreeType;
			
			/// Get the harmonic degree.
			HarmonicDegree harmonicDegree() const 
			{
				if (DerivedConstCast->empty()) 
				{
					return HarmonicDegree(0);
				}

				const typename Derived::const_iterator result(std::max_element(DerivedConstCast->begin(), DerivedConstCast->end(),
																               HarmonicDegreeBinaryPredicate<typename Derived::TermType>() ));

				return result->template get<HarmonicTermPosition>().harmonicDegree();
			}


			/// Get the harmonic order.
			/**
			 * The harmonic order is defined as the minimum harmonic degree of the terms composing the series.
			 */
			HarmonicDegree harmonicOrder() const 
			{
				if (DerivedConstCast->empty()) 
				{
					return HarmonicDegree(0);
				}

				const typename Derived::const_iterator result(std::min_element(DerivedConstCast->begin(), DerivedConstCast->end(),
							                                                   HarmonicOrderBinaryPredicate<typename Derived::TermType>() ));

				return result->template get<HarmonicTermPosition>().harmonicOrder();
			}


			/// Return true if all terms of the harmonic series are cosines.
			/**
			 * An empty series will return false.
			 */
			bool isCosine() const
			{
				typedef typename Derived::const_iterator const_iterator;

				if (DerivedConstCast->empty()) 
				{
					return false;
				}

				const const_iterator itf = DerivedConstCast->end();
				for (const_iterator it = DerivedConstCast->begin(); it != itf; ++it) 
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
			bool isSine() const
			{
				typedef typename Derived::const_iterator const_iterator;
				if (DerivedConstCast->empty()) 
				{
					return false;
				}
				const const_iterator itf = DerivedConstCast->end();
				for (const_iterator it = DerivedConstCast->begin(); it != itf; ++it) 
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
			HarmonicDegree basePartialHarmonicDegree(const PosTuple &posTuple) const 
			{
				if (DerivedConstCast->empty()) 
				{
					return HarmonicDegree(0);
				}
				const typename Derived::const_iterator result(std::max_element(DerivedConstCast->begin(), DerivedConstCast->end(),
							                                                   PartialHarmonicDegreeBinaryPredicate<typename Derived::TermType, PosTuple>(posTuple) ));
				
				return result->template get<HarmonicTermPosition>().partialHarmonicDegree(posTuple);
			}


			/// Get the mininum harmonic degree of the series for specific variables.
			template <class PosTuple>
			HarmonicDegree basePartialHarmonicOrder(const PosTuple &posTuple) const 
			{
				if (DerivedConstCast->empty()) 
				{
					return HarmonicDegree(0);
				}

				const typename Derived::const_iterator result(std::min_element(DerivedConstCast->begin(), DerivedConstCast->end(),
							                                                   PartialHarmonicOrderBinaryPredicate<typename Derived::TermType, PosTuple>(posTuple) ));
				
				return result->template get<HarmonicTermPosition>().partialHarmonicOrder(posTuple);
			}
	};
}

#undef DerivedConstCast
//#undef derived_cast

#endif
