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

#include "core/base_classes/expo_truncator.h"

namespace piranha
{
	// Static initialization for expo-based truncation.
	base_expo_truncator::container_type base_expo_truncator::m_expo_limits;

	void base_expo_truncator::clear_all()
	{
		m_expo_limits.clear();
	}

	void base_expo_truncator::print(std::ostream &stream)
	{
		if (m_expo_limits.empty()) {
			stream << "No exponent limits defined.";
		} else {
			const iterator it_f = m_expo_limits.end();
			for (iterator it = m_expo_limits.begin(); it != it_f;) {
				stream << it->first->name() << ',' << it->second;
				++it;
				if (it != it_f) {
					stream << '\n';
				}
			}
		}
	}

	std::string base_expo_truncator::__repr__()
	{
		std::ostringstream stream;
		print(stream);
		std::string retval(stream.str());
		return retval;
	}

	base_expo_truncator::iterator base_expo_truncator::find_argument(const psym_p &p)
	{
		const iterator it_f = m_expo_limits.end();
		iterator it(m_expo_limits.begin());
		for (; it != it_f; ++it) {
			if (it->first == p) {
				break;
			}
		}
		return it;
	}
}
