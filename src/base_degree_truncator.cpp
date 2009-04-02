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

#include <iostream>
#include <sstream>

#include "core/base_classes/degree_truncator.h"

namespace piranha
{
	// Static initialization for degree-based truncation.
	int base_degree_truncator::m_degree_limit = 0;
	bool base_degree_truncator::m_effective = false;

	void base_degree_truncator::unset()
	{
		m_effective = false;
	}

	void base_degree_truncator::print(std::ostream &stream)
	{
		if (m_effective) {
			stream << "Minimum degree limit: " << m_degree_limit;
		} else {
			stream << "No minimum degree limit set.";
		}
	}
}
