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

#ifndef PIRANHA_NORM_TRUNCATOR_H
#define PIRANHA_NORM_TRUNCATOR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <vector>

#include <boost/lambda/lambda.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "../config.h"
#include "../exceptions.h"
#include "../settings.h"
#include "../utils.h"

namespace piranha {
namespace truncators {

	/// Norm-based truncator.
	class PIRANHA_VISIBLE Norm
	{
			template <class ArgsTuple>
			class CompareNorm
			{
				public:
					CompareNorm(ArgsTuple const &argsTuple): argsTuple(argsTuple) {}

					template <class Term>
					bool operator()(Term const * const t1, Term const * const t2) const
					{
						double const norm1(t1->cf.norm(argsTuple) * t1->key.norm(argsTuple));
                        double const norm2(t2->cf.norm(argsTuple) * t2->key.norm(argsTuple));
						return norm1 > norm2;
					}

				private:
					const ArgsTuple &argsTuple;
			};


		public:
			template <class Series1, class Series2, class ArgsTuple>
			class GetType
			{
				public:
					typedef typename Series1::TermType TermType1;
					typedef typename Series2::TermType TermType2;
					typedef GetType Type;

					GetType(std::vector<TermType1 const *> &terms1, std::vector<TermType2 const *> &terms2, const ArgsTuple &argsTuple, bool initialise = true)
                        : terms1(terms1), terms2(terms2), argsTuple(argsTuple), term1(0), term2(0)
					{
						if (initialise) 
						{
							init();
						}
					}


					bool skip(TermType1 const **t1, TermType2 const **t2) const
					{
						return (term1 && t1 >= term1) || (term2 && t2 >= term2);
					}


					bool isEffective() const
					{
						return truncationLevel != 0.0;
					}


					// Returns the length of a development in powers of x that satisfies the condition that the
					// magnitude of the last term of the expansion with respect to x's magnitudes is m_truncation_level
					// times smaller.
					template <class Series, class ArgsTuple2>
					static std::size_t powerSeriesIterations(Series const &series,  int const start, int const stepSize,  ArgsTuple2 const &argsTuple)
					{
						// NOTE: share this check in some kind of base truncator class?
						if (stepSize < 1) 
						{
							PIRANHA_THROW(value_error, "Please use a step size of at least 1");
						}

						if (truncationLevel == 0.0) 
						{
							PIRANHA_THROW(value_error,  "No value set for norm-based truncation, "
								                        "cannot calculate limit of power series expansion");
						}

						const double norm = series.baseNorm(argsTuple);

						PIRANHA_ASSERT(norm >= 0.0);
						if (norm >= 1.0) 
						{
							PIRANHA_THROW(value_error, "The norm of the argument of the power series expansion is >= 1: "
								"the norm truncator is unable to give an estimate of the power series limit");
						}

						// Let's prevent log10(0) below.
						if (norm == 0.0) 
						{
							PIRANHA_THROW(value_error, "Unable to find a limit for the power series expansion of a series whose norm "
								"is zero");
						}

						int const retval = boost::numeric_cast<int>(std::ceil((std::log(truncationLevel) / std::log(norm) + stepSize - start) / stepSize + 1));

						// This could be negative if starting power is big enough. In this case return 0.
						return (retval >= 0) ? retval : 0;
					}



                    // return a vector of pointers to the series terms sorted according to the norm
					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::TermType const *> getSortedPointerVector(Series const &series, ArgsTuple2 const &argsTuple2)
					{
                        typedef std::vector<typename Series::TermType const *> RetValType; // a vector of pointers to series terms
						RetValType retval;

                        // create vector of pointers to the series terms in the retval vector
                        // the _1 means the first parameter of the functional operator (boost::lambda). The expresion determines the address (&) of the series term
                        // iterated through (series.begin()..series.end() and this address is inserted (insert_iterator) into the retval vector starting at retval.begin() 
						std::transform(series.begin(), series.end(), std::insert_iterator< RetValType >(retval, retval.begin()), &boost::lambda::_1);

						if (truncationLevel == 0) 
						{
							PIRANHA_THROW(value_error, "Cannot establish series ordering, norm truncator is not active");
						}

                        // sort according to the series term norm
						std::sort(retval.begin(), retval.end(), CompareNorm<ArgsTuple2>(argsTuple2));
						return retval;
					}

				protected:
					
					void init()
					{
						if (isEffective()) 
						{
							PIRANHA_ASSERT(terms1.size() >= 1 && terms2.size() >= 1);
							
                            // we need it twice. Construct only once
                            const CompareNorm<ArgsTuple> cmp(argsTuple);

							std::sort(terms1.begin(), terms1.end(), cmp);
							std::sort(terms2.begin(), terms2.end(), cmp);

							const double norm1 = calculateNorm(terms1);
                            const double norm2 = calculateNorm(terms2);
							// If one of the norms is zero, then we don't want to do any multiplication.
							const TermType1 **final1 = &(*(terms1.begin()));
							const TermType2 **final2 = &(*(terms2.begin()));

							if (norm1 == 0.0 || norm2 == 0.0) 
							{
								term1 = final1;
								term2 = final2;
								return;
							}

							double const limit1 = -norm2 + std::sqrt(norm2 * norm2 + truncationLevel * norm2 / norm1);
							double const limit2 = -norm1 + std::sqrt(norm1 * norm1 + truncationLevel * norm1 / norm2);
							double delta1 = 0.0;
							double delta2 = 0.0;
							// NOTE: the issue here is that it may happen that limit is greater than norm, in which case all the calculations
							// done for maximimzing the number of terms truncated do not apply. In such case for the time being we just turn
							// off truncation, but it should be probably dealt with differently.
							while (limit1 < norm1) 
							{
								const TermType1 **tmp = (term1) ? term1 - 1 : &(*(terms1.end() - 1));

								delta1 += (*tmp)->cf.norm(argsTuple) * (*tmp)->key.norm(argsTuple);
								
                                // If, by going to the next term, we exceed the delta1, leave term1 where it is and break out.
								if (delta1 >= limit1) 
								{
									break;
								}

								// Assign new limit.
								term1 = tmp;
								// If we reached the top, break out.
								if (tmp == final1) 
								{
									break;
								}
							}

							while (limit2 < norm2) 
							{
								const TermType2 **tmp = (term2) ? term2 - 1 : &(*(terms2.end() - 1));

								delta2 += (*tmp)->cf.norm(argsTuple) * (*tmp)->key.norm(argsTuple);
								
                                // If, by going to the next term, we exceed the delta2, leave term2 where it is and break out.
								if (delta2 >= limit2) 
								{
									break;
								}

								// Assign new limit.
								term2 = tmp;
								// If we reached the top, break out.
								if (tmp == final2) 
								{
									break;
								}
							}
						}
					}

				private:

					template <class Term>
					double calculateNorm(std::vector<Term const *> const &v) const
					{
						double retval = 0;
						std::size_t const end = v.size();
						for (std::size_t i = 0; i < end; ++i) 
						{
							retval += v[i]->cf.norm(argsTuple) * v[i]->key.norm(argsTuple);
						}
						return retval;
					}

					std::vector<TermType1 const *>	&terms1;
					std::vector<TermType2 const *>	&terms2;
					ArgsTuple const			        &argsTuple;
					TermType1 const		            **term1;
					TermType2 const		            **term2;
			};

			// Shared portion.
			static void set(double const &x)
			{
				if (x <= 0.0) 
				{
					PIRANHA_THROW(value_error, "Please use a positive number for norm truncation");
				}

				truncationLevel = x;
			}

			static void print(std::ostream &stream = std::cout);
			static void unset();

		private:

			static double truncationLevel;
	};
} }

#endif
