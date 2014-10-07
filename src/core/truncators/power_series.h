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
			class get_type:	public degree::template GetType<Series1, Series2, ArgsTuple>,
				            public norm::template   get_type<Series1, Series2, ArgsTuple>
			{
					typedef typename degree::template GetType<Series1, Series2, ArgsTuple> DegreeAncestor;
					typedef typename norm::template   get_type<Series1, Series2, ArgsTuple> NormAncestor;

					enum SelectedTruncator {degTruncator, normTruncator, nullTruncator};

				public:

					typedef get_type type;
					typedef typename Series1::term_type TermType1;
					typedef typename Series2::term_type TermType2;

					get_type(std::vector<TermType1 const *> &terms1, std::vector<TermType2 const *> &terms2, const ArgsTuple &argsTuple):
						DegreeAncestor(terms1, terms2, argsTuple, false), NormAncestor(terms1, terms2, argsTuple, false), activeTruncator(degTruncator)
					{
						DegreeAncestor::init();
					
                        if(!DegreeAncestor::isEffective()) 
						{
							NormAncestor::init();
						
                            if(!NormAncestor::isEffective()) 
							{
								activeTruncator = nullTruncator;

							} else 
							{
								activeTruncator = normTruncator;
							}

						} else 
						{
							activeTruncator = degTruncator;
						}
					}


					template <class T, class ArgsTuple2>
					static std::size_t powerSeriesIterations(const T &x, const int &start, const int &stepSize, const ArgsTuple2 &argsTuple) 
					{
						std::string msg("No useful truncation limit for a power series expansion could be "
							            "established by the power series truncator. The reported errors were:\n");

						try {
							return DegreeAncestor::powerSeriesIterations(x, start, stepSize, argsTuple);

						}catch (const value_error &ve) 
						{
							msg += std::string(ve.what()) + "\n";
						}

						try {
							return NormAncestor::powerSeriesIterations(x, start, stepSize, argsTuple);

						} catch (const value_error &ve) 
						{
							msg += std::string(ve.what()) + "\n";
						}

						piranha_throw(value_error,msg);
					}


					template <class Series, class ArgsTuple2>
					static std::vector<typename Series::term_type const *> getSortedPointerVector(const Series &s, const ArgsTuple2 &argsTuple)
					{
						std::string msg("The power series truncator was not able to establish a series ordering. The reported errors were:\n");
						try {
							return DegreeAncestor::getSortedPointerVector(s, argsTuple);

						} catch (const value_error &ve)
                        {
									msg += std::string(ve.what()) + "\n";
						}

						try {
							return NormAncestor::getSortedPointerVector(s, argsTuple);

						} catch (const value_error &ve)
                        {
									msg += std::string(ve.what()) + "\n";
						}

						piranha_throw(value_error,msg);
					}


					bool isEffective() const
                    {
						return activeTruncator != nullTruncator;
					}


					template <class T, class U>
					bool skip(const T &x1, const U &x2) const 
					{
						switch (activeTruncator)
                        {
							case degTruncator:	 return DegreeAncestor::skip(x1,x2);
						
                            case normTruncator:  return NormAncestor::skip(x1,x2);
							
                            case nullTruncator:  piranha_assert(false);
						}

						piranha_assert(false);
						return false;
					}

				private:
					
					SelectedTruncator activeTruncator;
			};
	};
}
}

#endif
