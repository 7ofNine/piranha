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

#ifndef PIRANHA_DEGREE_TRUNCATOR_H
#define PIRANHA_DEGREE_TRUNCATOR_H

#include <algorithm> // For sorting.
#include <cmath> // For std::ceil.
#include <cstddef>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

#include <boost/lambda/lambda.hpp>
#include <boost/tuple/tuple.hpp>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "../settings.h" // For debug messages.

namespace piranha { 
namespace truncators {

	/// Truncator based on the minium degree of the series.
	class PIRANHA_VISIBLE Degree
	{

			enum TruncationMode {
				TRUNCATION_DEGREE         = 0,
				TRUNCATION_PARTIAL_DEGREE = 1,
				TRUNCATION_INACTIVE       = 2,
			};


            // comparison functor for  Exponent Vectors
			template <int ExpoTermPos>
			class CompareOrder
			{
                public: 

				template <class Term>
				bool operator()(Term const * const t1, Term const * const t2) const
				{
					typedef typename Term::template Component<ExpoTermPos>::Type::DegreeType DegreeType;

					DegreeType const md1(t1->template get<ExpoTermPos>().order());
                    DegreeType const md2(t2->template get<ExpoTermPos>().order());

					return md1 < md2;
				}
			};


            //comparison functor for partial i.e. subset of arguments Exponent Vector
            //the position tuple picks the element in the term keys. 
            // The position has to be determind beofre this can be used
			template <int ExpoTermPos, class PositionTuple>
			class ComparePartialOrder
			{
                public:

				ComparePartialOrder(const PositionTuple &positionTuple): positionTuple(positionTuple) {}

				template <class Term>
				bool operator()(Term const * const t1, Term const * const t2) const
				{
					typedef typename Term::template Component<ExpoTermPos>::Type::DegreeType DegreeType;

					DegreeType const md1(t1->template get<ExpoTermPos>().partialOrder(positionTuple));
					DegreeType const md2(t2->template get<ExpoTermPos>().partialOrder(positionTuple));

                    return md1 < md2;
				}

				PositionTuple const &positionTuple;
			};


		public:

			template <class Series1, class Series2, class ArgsTuple>
			class GetType
			{
					typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type PosTupleType;

					static const int exponentTermPosition = Series1::exponentTermPosition;
					static const int exponentArgsPosition = Series1::exponentArgsPosition;

				public:
					typedef typename Series1::TermType TermType1;
					typedef typename Series2::TermType TermType2;
					typedef GetType type;

					GetType(std::vector<TermType1 const *> &terms1, std::vector<TermType2 const *> &terms2, const ArgsTuple &argsTuple, bool initialise = true)
                            :terms1(terms1), terms2(terms2), argsTuple(argsTuple)
					{
						// Some static checks.
                        static_assert(Series1::exponentArgsPosition == Series2::exponentArgsPosition, "");
                        static_assert(Series1::exponentTermPosition == Series2::exponentTermPosition, "");

						// Convert psyms vector into position tuple only if we are truncating to partial degree.
						if (truncationMode == TRUNCATION_PARTIAL_DEGREE) 
						{
							positionTuple = psyms2pos(psyms, argsTuple);
						}

						if (initialise) 
						{
							init();
						}
					}


					bool skip(TermType1 const **t1, TermType2 const **t2) const
					{
						switch (truncationMode)
                        {
							case TRUNCATION_DEGREE:         return   ((*t1)->template get<exponentTermPosition>().order() + (*t2)->template get<exponentTermPosition>().order() >= degreeLimit);

							case TRUNCATION_PARTIAL_DEGREE: return   ((*t1)->template get<exponentTermPosition>().partialOrder(positionTuple) 
                                                                    + (*t2)->template get<exponentTermPosition>().partialOrder(positionTuple) >= degreeLimit);

							case TRUNCATION_INACTIVE:       PIRANHA_ASSERT(false); // We should never get there.
				 	    }

						return false;
					}


					// Number of a iterations of a power series development of a power series.
					// NOTE: if start is negative, it is assumed that negative powers of the input series
					// have a minimum degree which is proportional to the input series' and with its sign changed.
					template <class PowerSeries, class ArgsTuple2>
					static std::size_t powerSeriesIterations(PowerSeries const &s, int const start, int const stepSize, ArgsTuple2 const &argsTuple2)
					{
						if (stepSize < 1) 
						{
							PIRANHA_THROW(value_error, "Please use a step size of at least 1");
						}

						if (truncationMode == TRUNCATION_INACTIVE) 
						{
							PIRANHA_THROW(value_error, "Cannot calculate the limit of a power series expansion "
								                       "if no degree limit has been set");
						}

						if (s.empty()) 
						{
							PIRANHA_THROW(value_error, "Cannot calculate the limit of the power series expansion of "
								                       "an empty power series");
						}

						// order will be either total or partial, depending on the mode.
						mp_rational order(0);

						switch (truncationMode)
                        {
							case TRUNCATION_DEGREE:         order = s.order();
								                            break;

							case TRUNCATION_PARTIAL_DEGREE: order = s.basePartialOrder(psyms2pos(psyms, argsTuple2));
								                            break;

							case TRUNCATION_INACTIVE:       PIRANHA_ASSERT(false);  // should not happen. An error
						}

						if (order <= 0) 
						{
							PIRANHA_THROW(value_error, "Can not calculate the limit of a power series expansion if the (partial) minimum degree of the series is negative or zero");
						}

						if (degreeLimit < 0) 
						{
							PIRANHA_THROW(value_error, "Can not calculate the limit of a power series expansion if the minimum degree limit is negative");
						}


						// (mp_rational(limit) / order - start) / step_size + 1;
						mp_rational tmp(degreeLimit);
						tmp /= order;
						tmp -= start;
						tmp /= stepSize;
						tmp += 1;

						if (tmp > 0) 
						{
							// Take the floor of tmp and convert to integer.
							const int retval = (tmp.get_num() / tmp.get_den()).to_int();

							if (tmp == retval) 
							{
								// If tmp was an integer in the beginning, we want to take the number
								// immediately preceding it (or zero).
								return (retval > 0) ? (retval - 1) : 0;

							} else 
							{
								// If tmp was not an integer, let's take the floor.
								return retval;
							}
						} else 
						{
							PIRANHA_DEBUG(std::cout << "Negative power series limit calculated, inserting 0 instead." << '\n');
							return 0;
						}
					}


                    // return a vector of pointers to the series terms sorted according to the degree
					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::TermType const *> getSortedPointerVector(const Series &s, const ArgsTuple2 &argsTuple2)
					{
						
                        typedef std::vector<typename Series::TermType const *> RetValType;
                        RetValType retval;
                        
                        // create vector of pointers to the series terms in the retval vector
                        // the _1 means the first parameter of the functional operator (boost::lambda). The expresion determines the address (&) of the series term
                        // iterated through (series.begin()..series.end() and this address is inserted (insert_iterator) into the retval vector starting at retval.begin() 
                        std::transform(s.begin(), s.end(), std::insert_iterator< RetValType >(retval, retval.begin()),	&boost::lambda::_1);

						switch (truncationMode)
                        {
							case TRUNCATION_DEGREE:         std::sort(retval.begin(), retval.end(), CompareOrder<Series::exponentTermPosition>());
								                            break;

							case TRUNCATION_PARTIAL_DEGREE:	{
								                                typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple2>::value>::Type PosTupleType;
								                                PosTupleType const pos_tuple(psyms2pos(psyms, argsTuple2));

								                                if (pos_tuple.template get<Series::exponentArgsPosition>().size() > 0) 
								                                {
									                                std::sort(retval.begin(), retval.end(), ComparePartialOrder<Series::exponentTermPosition, PosTupleType>(pos_tuple));

								                                } else
								                                {
									                                PIRANHA_THROW(value_error, "Cannot establish series ordering, partial degree truncator is not effective on this series");
								                                }
								                            }
                                                            break;

							case TRUNCATION_INACTIVE:       PIRANHA_THROW(value_error, "cannot establish series ordering, degree truncator is not active");
						}

						return retval;
					}


                    // test if we are truncating
					bool isEffective() const
					{
						switch (truncationMode)
                        {
							case TRUNCATION_INACTIVE:       return false;

							case TRUNCATION_DEGREE:         return true;
							
                            case TRUNCATION_PARTIAL_DEGREE: return (positionTuple.template get<exponentArgsPosition>().size() > 0);

								           // In case of partial degree truncator is effective only if the position tuple
								           // contains some elements.
								
						}

						PIRANHA_ASSERT(false);
						return false;  // shut up compiler. will never et here
					}

				protected:

					void init()
					{
						switch (truncationMode)
                        {
							case TRUNCATION_DEGREE:         std::sort(terms1.begin(), terms1.end(), CompareOrder<exponentTermPosition>());
								                            std::sort(terms2.begin(), terms2.end(), CompareOrder<exponentTermPosition>());
								                            break;

							case TRUNCATION_PARTIAL_DEGREE:	// We need to do the sorting only if the position tuple
								                            // contains some elements.
								                            if (positionTuple.template get<exponentArgsPosition>().size() > 0) 
								                            {
									                            std::sort(terms1.begin(), terms1.end(), ComparePartialOrder<exponentTermPosition, PosTupleType>(positionTuple));
                                                                std::sort(terms2.begin(), terms2.end(), ComparePartialOrder<exponentTermPosition, PosTupleType>(positionTuple));
								                            }
								                            break;
							case TRUNCATION_INACTIVE: ;   // nothing todo

						}
					}

				private:

					std::vector<TermType1 const *>	&terms1;
					std::vector<TermType2 const *>	&terms2;
					const ArgsTuple			        &argsTuple;
					PosTupleType			         positionTuple;
			};


			static void set(const int);
			static void set(const mp_rational &);
			static void set(const std::string &, const int);
			static void set(const std::string &, const mp_rational &);
			static void set(const std::vector<std::string> &, const int);
			static void set(const std::vector<std::string> &, const mp_rational &);
			static void unset();
			static void print(std::ostream &stream = std::cout);

		private:

			static mp_rational	    degreeLimit;
			static VectorPsym	    psyms;
			static TruncationMode   truncationMode;
	};
} }

#endif
