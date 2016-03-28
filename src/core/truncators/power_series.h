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

#ifndef PIRANHA_POWER_SERIES_TRUNCATOR
#define PIRANHA_POWER_SERIES_TRUNCATOR

#include <cstddef>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "degree.h"
#include "norm.h"

namespace piranha
{
namespace truncators
{
	/// Truncator for power series.
	class PowerSeries
	{
		public:

			template <class Series1, class Series2, class ArgsTuple>
			class GetType:	public Degree::template GetType<Series1, Series2, ArgsTuple>,
				            public Norm::template   GetType<Series1, Series2, ArgsTuple>
			{
					typedef typename Degree::template GetType<Series1, Series2, ArgsTuple> DegreeAncestor;
					typedef typename Norm::template   GetType<Series1, Series2, ArgsTuple> NormAncestor;

					enum SelectedTruncator {
                        TRUNCATOR_NULL   = 0,
                        TRUNCATOR_DEGREE = 1,
                        TRUNCATOR_NORM   = 2
                        };

				public:

					typedef GetType Type;
					typedef typename Series1::TermType TermType1;
					typedef typename Series2::TermType TermType2;

					GetType(std::vector<TermType1 const *> &terms1, std::vector<TermType2 const *> &terms2, const ArgsTuple &argsTuple)
                        : DegreeAncestor(terms1, terms2, argsTuple, false), NormAncestor(terms1, terms2, argsTuple, false), activeTruncator(TRUNCATOR_DEGREE)
					{
						DegreeAncestor::init();
					
                        if(!DegreeAncestor::isEffective()) 
						{
							NormAncestor::init();
						
                            if(!NormAncestor::isEffective()) 
							{
								activeTruncator = TRUNCATOR_NULL;

							} else 
							{
								activeTruncator = TRUNCATOR_NORM;
							}

						} else 
						{
							activeTruncator = TRUNCATOR_DEGREE;
						}
					}


					template <class Series, class ArgsTuple2>
					static std::size_t powerSeriesIterations(const Series &series, const int &start, const int &stepSize, const ArgsTuple2 &argsTuple2) 
					{
						std::string msg("No useful truncation limit for a power series expansion could be "
							            "established by the power series truncator. The reported errors were:\n");

						try {
							return DegreeAncestor::powerSeriesIterations(series, start, stepSize, argsTuple2);

						}catch (const value_error &ve) 
						{
							msg += std::string(ve.what()) + "\n";
						}

						try {
							return NormAncestor::powerSeriesIterations(series, start, stepSize, argsTuple2);

						} catch (const value_error &ve) 
						{
							msg += std::string(ve.what()) + "\n";
						}

						PIRANHA_THROW(value_error, msg);
					}


					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::TermType const *> getSortedPointerVector(const Series &series, const ArgsTuple2 &argsTuple2)
					{
						std::string msg("The power series truncator was not able to establish a series ordering. The reported errors were:\n");
						try {
							return DegreeAncestor::getSortedPointerVector(series, argsTuple2);

						} catch (const value_error &ve)
                        {
									msg += std::string(ve.what()) + "\n";
						}

						try {
							return NormAncestor::getSortedPointerVector(series, argsTuple2);

						} catch (const value_error &ve)
                        {
									msg += std::string(ve.what()) + "\n";
						}

						PIRANHA_THROW(value_error,msg);
					}


					bool isEffective() const
                    {
						return activeTruncator != TRUNCATOR_NULL;
					}


					template <class TermType1, class TermType2>
					bool skip(const TermType1 &x1, const TermType2 &x2) const 
					{
						switch (activeTruncator)
                        {
							case TRUNCATOR_DEGREE:  return DegreeAncestor::skip(x1, x2);
						
                            case TRUNCATOR_NORM:    return NormAncestor::skip(x1, x2);
							
                            case TRUNCATOR_NULL:    PIRANHA_ASSERT(false);
						}

						PIRANHA_ASSERT(false);
						return false;
					}

				private:
					
					SelectedTruncator activeTruncator;
			};
	};
}
}

#endif
