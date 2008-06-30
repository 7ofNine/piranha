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
#include <string>

#include "core/base_classes/norm_truncator.h"
#include "core/stream_manager.h"

namespace piranha
{
// Initial value for norm-based truncation.
	double base_norm_truncator::m_truncation_level = 1E-6;

	void base_norm_truncator::print(std::ostream &stream)
	{
		stream_manager::setup_print(stream);
		stream << "Truncation level: " << m_truncation_level;
	}

	std::string base_norm_truncator::py_repr()
	{
		std::ostringstream stream;
		print(stream);
		std::string retval(stream.str());
		return retval;
	}
}
