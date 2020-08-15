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

#ifndef PIRANHA_BASE_POWER_SERIES_H
#define PIRANHA_BASE_POWER_SERIES_H

#include <algorithm> // For max_element and min_element.

#include "../config.h"
#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	/// Power series toolbox.
	template <int ExpoArgsPosition, int ExpoTermPosition, class Degree, class Derived>
	class BasePowerSeries
	{
        static_assert(ExpoArgsPosition >= 0, "Invalid exponent argument position.");


            //compare functor for the total exponential degree
			template <class Term>
			class CompareDegree
			{
                public:
				explicit CompareDegree(VectorPsym const & s) :symbols(s) {}

				bool operator()(Term const &term1, Term const &term2) const 
                {
					return (term1.template get<ExpoTermPosition>().degree(symbols) < term2.template get<ExpoTermPosition>().degree(symbols));
				}

				private:
					VectorPsym const & symbols;
			};


            // comparison functor dor the exponential order 
			template <class Term>
			class CompareOrder
			{
				public:
					explicit CompareOrder(VectorPsym const & s) :symbols(s) {}

					bool operator()(Term const &term1, Term const &term2) const 
					{
						return (term1.template get<ExpoTermPosition>().order(symbols) < term2.template get<ExpoTermPosition>().order(symbols));
					}

				private:
					VectorPsym const & symbols;
			};


			template <class Term>
			class XCompareOrder      //GUT possible fix for parameter free order
			{
			public:
				explicit XCompareOrder(){}

				bool operator()(Term const& term1, Term const& term2) const
				{
					return (term1.template get<ExpoTermPosition>().xorder() < term2.template get<ExpoTermPosition>().xorder());
				}

			};


            // comparison functor for partial degree. Only arguments that are marked by the positionTuple are considered
			template <class Term, class PositionTuple>
			class  ComparePartialDegree
			{
                public:

				ComparePartialDegree(PositionTuple const &p):positionTuple(p) {}

				bool operator()(Term const &term1, Term const &term2) const 
                {
					return (term1.template get<ExpoTermPosition>().partialDegree(positionTuple) < term2.template get<ExpoTermPosition>().partialDegree(positionTuple));
				}

                private:

				const PositionTuple &positionTuple;
			};


			template <class Term, class PositionTuple>
			class ComparePartialOrder
			{
                public:

				ComparePartialOrder(PositionTuple const &p):positionTuple(p) {}

				bool operator()(Term const &term1, Term const &term2) const 
                {
					return (term1.template get<ExpoTermPosition>().partialOrder(positionTuple) < term2.template get<ExpoTermPosition>().partialOrder(positionTuple));
				}

                private:

				const PositionTuple &positionTuple;
			};


		public:

			static const int exponentArgsPosition = ExpoArgsPosition;
			static const int exponentTermPosition = ExpoTermPosition;
			typedef Degree DegreeType;

			/// Get the degree of the power series. degree is the maximum degree of the exponent vector
			Degree degree(VectorPsym const & symbols) const
            {
				if (derived_const_cast->empty()) 
                {
					return Degree(0);
				}

				const typename Derived::const_iterator result(std::max_element(derived_const_cast->begin(), derived_const_cast->end(), CompareDegree<typename Derived::TermType>(symbols) ));

				return result->template get<ExpoTermPosition>().degree(symbols);
			}


			/// Get the order of the power series.
			/**
			 * The order is defined as the minimum degree of the terms composing the series.
			 */
			Degree order(VectorPsym const & symbols) const 
            {
				if (derived_const_cast->empty()) 
                {
					return Degree(0);
				}

				const typename Derived::const_iterator result(std::min_element(derived_const_cast->begin(), derived_const_cast->end(), CompareOrder<typename Derived::TermType>(symbols) ));
				return result->template get<ExpoTermPosition>().order(symbols);
			}

			Degree xorder() const
			{
				if (derived_const_cast->empty())
				{
					return Degree(0);
				}
				const typename Derived::const_iterator result(std::min_element(derived_const_cast->begin(), derived_const_cast->end(), XCompareOrder<typename Derived::TermType>()));
				return result->template get<ExpoTermPosition>().xorder();


			}

		//protected:
			/// Get the degree of the power series for specific variables.
			template <class PositionTuple>
			Degree basePartialDegree(PositionTuple const &positionTuple) const 
            {
				if (derived_const_cast->empty()) 
                {
					return Degree(0);
				}

				auto begin = derived_const_cast->begin();
				auto end = derived_const_cast->end();

				//for (Derived::const_iterator it = begin; it ++; it != end)
				//{
				//	auto it2 = it;
				//}
				//const typename Derived::const_iterator result(std::max_element(derived_const_cast->begin(), derived_const_cast->end(),
					const typename Derived::const_iterator result(std::max_element(begin, end, 
                                                              ComparePartialDegree<typename Derived::TermType, PositionTuple>(positionTuple) ));

				return result->template get<ExpoTermPosition>().partialDegree(positionTuple);
			}


			/// Get the mininum degree of the power series for specific variables.
			template <class PositionTuple>
			Degree basePartialOrder(PositionTuple const &positionTuple) const 
            {
				if (derived_const_cast->empty()) 
                {
					return Degree(0);
				}

				const typename Derived::const_iterator result(std::min_element(derived_const_cast->begin(), derived_const_cast->end(),
							                                  ComparePartialOrder<typename Derived::TermType, PositionTuple>(positionTuple) ));

				return result->template get<ExpoTermPosition>().partialOrder(positionTuple);
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
