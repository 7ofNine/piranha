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

#include <boost/lexical_cast.hpp>
#include <iostream>
#include <sstream>
#include <string>

#include "core/exceptions.h"
#include "core/settings.h"
#include "core/truncators/norm.h"

using namespace piranha::truncators;

namespace piranha
{
    // allocation for static member and
	// and initial value for norm-based truncation.
    //
	double Norm::truncationLevel = 0.0;

    // print the current status of the norm truncator
	void Norm::print(std::ostream &stream)
	{
		if (truncationLevel > 0)
         {
			stream << "Truncation level: " << boost::lexical_cast<std::string>(truncationLevel);
		} else
        {
			stream << "No truncation level set.";
		}
	}


	void Norm::unset()
	{
		truncationLevel = 0.0;
	}
}
