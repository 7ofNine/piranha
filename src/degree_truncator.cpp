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
	// Static initialization for degree-based truncation.
	mp_rational            degree::degreeLimit;
	degree::TruncationMode degree::truncationMode = degree::TruncationInactive;
	VectorPsym             degree::psyms;

	void degree::unset()
	{
		truncationMode = TruncationInactive;
	}

	void degree::set(const int n)
	{
		degreeLimit    = n;
		truncationMode = TruncationDeg;
	}

	void degree::set(const mp_rational &r)
	{
		degreeLimit    = r;
		truncationMode = TruncationDeg;
	}

	void degree::set(const std::string &s, const int n)
	{
		degreeLimit    = n;
		psyms          = names2psyms(std::vector<std::string>(1, s)); // setup input string as vector of strings
		truncationMode = TruncationPartialDeg;
	}

	void degree::set(const std::string &s, const mp_rational &r)
	{
		degreeLimit    = r;
		psyms          = names2psyms(std::vector<std::string>(1, s)); // setup input string as vector of strings
		truncationMode = TruncationPartialDeg;
	}

	void degree::set(const std::vector<std::string> &vs, const int n)
	{
		if (!vs.size())
        {
			set(n);
			return;
		}
		degreeLimit    = n;
		psyms          = names2psyms(vs);
		truncationMode = TruncationPartialDeg;
	}

	void degree::set(const std::vector<std::string> &vs, const mp_rational &r)
	{
		if (!vs.size())
        {
			set(r);
			return;
		}
		degreeLimit     = r;
		psyms           = names2psyms(vs);
		truncationMode  = TruncationPartialDeg;
	}

void degree::print(std::ostream &stream)
{
	switch (truncationMode)
    {
		case TruncationInactive: stream << "No degree limit set.";
			                     break;

		case TruncationDeg:      stream << "Degree limit: " << degreeLimit;
			                     break;

		case TruncationPartialDeg:  {   stream << "Partial degree limit: " << degreeLimit << ", Affected symbols: [";
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

