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

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "core/mp.h"
#include "core/psym.h"
#include "core/truncators/degree.h"

using namespace piranha::truncators;

namespace piranha
{ 
	// Static initialization and storage allocation for degree-based truncation.
	mp_rational            Degree::degreeLimit;
	Degree::TruncationMode Degree::truncationMode = Degree::TRUNCATION_INACTIVE;
	VectorPsym             Degree::psyms;

	void Degree::unset()
	{
		truncationMode = TRUNCATION_INACTIVE;
	}

	void Degree::set(int const n)
	{
		degreeLimit    = n;
		truncationMode = TRUNCATION_DEGREE;
	}

	void Degree::set(mp_rational const &r)
	{
		degreeLimit    = r;
		truncationMode = TRUNCATION_DEGREE;
	}


    // set partial truncation for name s at degree n
	void Degree::set(std::string const &s, const int n)
	{
		degreeLimit    = n;
		psyms          = names2psyms(std::vector<std::string>(1, s)); // setup input string as vector of strings
		truncationMode = TRUNCATION_PARTIAL_DEGREE;
	}

    // set partial truncation for name s at rational degree r
	void Degree::set(std::string const &s, mp_rational const &r)
	{
		degreeLimit    = r;
		psyms          = names2psyms(std::vector<std::string>(1, s)); // setup input string as vector of strings
		truncationMode = TRUNCATION_PARTIAL_DEGREE;
	}

    // set partial truncation for several names in vs at degree n.
    // they all have the same degree
	void Degree::set(std::vector<std::string> const &vs, int const n)
	{
		if (!vs.size())
        {
			set(n);
			return;
		}

		degreeLimit    = n;
		psyms          = names2psyms(vs);
		truncationMode = TRUNCATION_PARTIAL_DEGREE;
	}

    // set partial truncation for several names in vs at rational degree r.
    // they all have the same degree
	void Degree::set(const std::vector<std::string> &vs, const mp_rational &r)
	{
		if (!vs.size())
        {
			set(r);
			return;
		}

		degreeLimit     = r;
		psyms           = names2psyms(vs);
		truncationMode  = TRUNCATION_PARTIAL_DEGREE;
	}

    // print current degree truncation settings
    void Degree::print(std::ostream &stream)
    {
	    switch (truncationMode)
        {
		    case TRUNCATION_INACTIVE:        stream << "No degree limit set.";
			                                 break;

		    case TRUNCATION_DEGREE:          stream << "Degree limit: " << degreeLimit;
			                                 break;

		    case TRUNCATION_PARTIAL_DEGREE:  {
                                                stream << "Partial degree limit: " << degreeLimit << ", Affected symbols: [";
			                                    for (std::size_t i = 0; i < psyms.size(); ++i)
                                                {
				                                    stream << '\'' << psyms[i].getName() << '\'';

                                                    if (i < psyms.size() - 1)
                                                    {
					                                    stream << " ";
				                                    }
			                                    }
			                                    stream << ']';
	                                         }
        }
    }
}

