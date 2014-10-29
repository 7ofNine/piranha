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
#include <boost/lambda/lambda.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath> // For std::ceil.
#include <cstddef>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "../settings.h" // For debug messages.

namespace piranha { 
namespace truncators {

	/// Truncator based on the minium degree of the series.
	class __PIRANHA_VISIBLE degree
	{

			enum TruncationMode {
				TruncationDeg        = 0,
				TruncationPartialDeg = 1,
				TruncationInactive   = 2,
			};


			template <int ExpoTermPos>
			struct order_comparison
			{
				template <class Term>
				bool operator()(const Term *t1, const Term *t2) const
				{
					typedef typename Term::template Component<ExpoTermPos>::Type::degree_type DegreeType;

					DegreeType const md1(t1->template get<ExpoTermPos>().order());
                    DegreeType const md2(t2->template get<ExpoTermPos>().order());

					return md1 < md2;
				}
			};


			template <int ExpoTermPos, class PosTuple>
			struct partial_order_comparison
			{
				partial_order_comparison(const PosTuple &posTuple): posTuple(posTuple) {}

				template <class Term>
				bool operator()(const Term *t1, const Term *t2) const
				{
					typedef typename Term::template Component<ExpoTermPos>::Type::degree_type DegreeType;

					DegreeType const md1(t1->template get<ExpoTermPos>().partialOrder(posTuple));
					DegreeType const md2(t2->template get<ExpoTermPos>().partialOrder(posTuple));

                    return md1 < md2;
				}

				PosTuple const &posTuple;
			};


		public:

			template <class Series1, class Series2, class ArgsTuple>
			class GetType
			{
					typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type PosTupleType;

					static const int expo_term_pos = Series1::expo_term_position;
					static const int expo_args_pos = Series1::expo_args_position;

				public:
					typedef typename Series1::TermType TermType1;
					typedef typename Series2::TermType TermType2;
					typedef GetType type;

					GetType(std::vector<TermType1 const *> &t1, std::vector<TermType2 const *> &t2, const ArgsTuple &argsTuple, bool initialise = true)
                            :m_t1(t1), m_t2(t2), m_argsTuple(argsTuple)
					{
						// Some static checks.
						PIRANHA_STATIC_CHECK(Series1::expo_args_position == Series2::expo_args_position, "");
						PIRANHA_STATIC_CHECK(Series1::expo_term_position == Series2::expo_term_position, "");

						// Convert psyms vector into position tuple only if we are truncating to partial degree.
						if (truncationMode == TruncationPartialDeg) 
						{
							m_pos_tuple = psyms2pos(psyms, m_argsTuple);
						}
						if (initialise) 
						{
							init();
						}
					}


					bool skip(const TermType1 **t1, const TermType2 **t2) const
					{
						switch (truncationMode)
                        {
							case TruncationDeg:      return ((*t1)->template get<expo_term_pos>().order() + (*t2)->template get<expo_term_pos>().order() >= degreeLimit);

							case TruncationPartialDeg:	   return ((*t1)->template get<expo_term_pos>().partialOrder(m_pos_tuple) + (*t2)->template get<expo_term_pos>().partialOrder(m_pos_tuple) >= degreeLimit);

							case TruncationInactive: PIRANHA_ASSERT(false); // We should never get there.
				 	    }
						return false;
					}


					// Number of a iterations of a power series development of a power series.
					// NOTE: if start is negative, it is assumed that negative powers of the input series
					// have a minimum degree which is proportional to the input series' and with its sign changed.
					template <class PowerSeries, class ArgsTuple2>
					static std::size_t powerSeriesIterations(const PowerSeries &s, const int &start, const int &stepSize, const ArgsTuple2 &argsTuple)
					{
						if (stepSize < 1) 
						{
							PIRANHA_THROW(value_error,"please use a step size of at least 1");
						}

						if (truncationMode == TruncationInactive) 
						{
							PIRANHA_THROW(value_error,"cannot calculate the limit of a power series expansion "
								"if no degree limit has been set");
						}

						if (s.empty()) 
						{
							PIRANHA_THROW(value_error,"cannot calculate the limit of the power series expansion of "
								"an empty power series");
						}

						// order will be either total or partial, depending on the mode.
						mp_rational order(0);

						switch (truncationMode)
                        {
							case TruncationDeg:   order = s.order();
								                  break;

							case TruncationPartialDeg: order = s.base_partial_order(psyms2pos(psyms, argsTuple));
								        break;

							case TruncationInactive: PIRANHA_ASSERT(false);
						}

						if (order <= 0) 
						{
							PIRANHA_THROW(value_error, "cannot calculate the limit of a power series expansion if the (partial) minimum degree of the series is negative or zero");
						}

						if (degreeLimit < 0) 
						{
							PIRANHA_THROW(value_error, "cannot calculate the limit of a power series expansion if the minimum degree limit is negative");
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
							__PDEBUG(std::cout << "Negative power series limit calculated, inserting 0 instead." << '\n');
							return 0;
						}
					}


					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::TermType const *> getSortedPointerVector(const Series &s, const ArgsTuple2 &argsTuple)
					{
						std::vector<typename Series::TermType const *> retval;
						std::transform(s.begin(), s.end(), std::insert_iterator<std::vector<typename Series::TermType const *> >(retval, retval.begin()),	&boost::lambda::_1);

						switch (truncationMode)
                        {
							case TruncationDeg:   std::sort(retval.begin(), retval.end(), order_comparison<Series::expo_term_position>());
								        break;

							case TruncationPartialDeg:	{
								            typedef typename NTuple<std::vector<std::pair<bool,std::size_t> >, boost::tuples::length<ArgsTuple2>::value>::Type PosTupleType;
								            PosTupleType const pos_tuple(psyms2pos(psyms, argsTuple));

								            if (pos_tuple.template get<Series::expo_args_position>().size() > 0) 
								            {
									            std::sort(retval.begin(), retval.end(),partial_order_comparison<Series::expo_term_position,PosTupleType>(pos_tuple));

								            } else
								            {
									            PIRANHA_THROW(value_error,"cannot establish series ordering, partial degree truncator is not effective on this series");
								            }
								        }
                                        break;

							case TruncationInactive: PIRANHA_THROW(value_error,"cannot establish series ordering, degree truncator is not active");
						}

						return retval;
					}


					bool isEffective() const
					{
						switch (truncationMode)
                        {
							case TruncationInactive: return false;

							case TruncationDeg:      return true;
							
                            case TruncationPartialDeg:    return (m_pos_tuple.template get<expo_args_pos>().size() > 0);

								           // In case of partial degree truncator is effective only if the position tuple
								           // contains some elements.
								
						}

						PIRANHA_ASSERT(false);
						return false;
					}

				protected:

					void init()
					{
						switch (truncationMode)
                        {
							case TruncationDeg:   std::sort(m_t1.begin(), m_t1.end(), order_comparison<expo_term_pos>());
								        std::sort(m_t2.begin(), m_t2.end(), order_comparison<expo_term_pos>());
								        break;

							case TruncationPartialDeg:	// We need to do the sorting only if the position tuple
								        // contains some elements.
								        if (m_pos_tuple.template get<expo_args_pos>().size() > 0) 
								        {
									        std::sort(m_t1.begin(), m_t1.end(), partial_order_comparison<expo_term_pos, PosTupleType>(m_pos_tuple));
                                            std::sort(m_t2.begin(), m_t2.end(), partial_order_comparison<expo_term_pos, PosTupleType>(m_pos_tuple));
								        }
								        break;
							case TruncationInactive: ;

						}
					}

				private:

					std::vector<TermType1 const *>	&m_t1;
					std::vector<TermType2 const *>	&m_t2;
					const ArgsTuple			        &m_argsTuple;
					PosTupleType			         m_pos_tuple;
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
