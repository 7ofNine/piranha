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

#include "core/base_classes/degree_truncator.h"
#include "core/psym.h"

namespace piranha
{
	// Static initialization for degree-based truncation.
	int degree_truncator::m_degree_limit = 0;
	degree_truncator::mode degree_truncator::m_mode = degree_truncator::inactive;
	vector_psym degree_truncator::m_psyms;

	void degree_truncator::unset()
	{
		m_mode = inactive;
	}

	void degree_truncator::set(const int &n) {
		m_degree_limit = n;
		m_mode = deg;
	}

	void degree_truncator::set(const vector_psym &v, const int &n) {
		m_degree_limit = n;
		m_psyms = v;
		m_mode = p_deg;
	}

	void degree_truncator::print(std::ostream &stream)
	{
		switch (m_mode) {
			case inactive:
				stream << "No minimum degree limit set.";
				break;
			case deg:
				stream << "Minimum degree limit: " << m_degree_limit;
				break;
			case p_deg:
				stream << "Partial minimum degree limit: " << m_degree_limit << '\n';
				stream << "Affected symbols: [";
				for (size_t i = 0; i < m_psyms.size(); ++i) {
					stream << m_psyms[i].get_name();
					if (i < m_psyms.size() - 1) {
						stream << " ";
					}
				}
				stream << ']';
		}
	}
}
