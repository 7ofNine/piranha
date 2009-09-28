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
	mp_rational degree::m_degree_limit;
	degree::mode degree::m_mode = degree::inactive;
	vector_psym degree::m_psyms;

	void degree::unset()
	{
		m_mode = inactive;
	}

	void degree::set(const int &n)
	{
		m_degree_limit = n;
		m_mode = deg;
	}

	void degree::set(const mp_rational &r)
	{
		m_degree_limit = r;
		m_mode = deg;
	}

	void degree::set(const std::string &s, const int &n)
	{
		m_degree_limit = n;
		m_psyms = names2psyms(std::vector<std::string>(1,s));
		m_mode = p_deg;
	}

	void degree::set(const std::string &s, const mp_rational &r)
	{
		m_degree_limit = r;
		m_psyms = names2psyms(std::vector<std::string>(1,s));
		m_mode = p_deg;
	}

	void degree::set(const std::vector<std::string> &vs, const int &n)
	{
		if (!vs.size()) {
			set(n);
			return;
		}
		m_degree_limit = n;
		m_psyms = names2psyms(vs);
		m_mode = p_deg;
	}

	void degree::set(const std::vector<std::string> &vs, const mp_rational &r)
	{
		if (!vs.size()) {
			set(r);
			return;
		}
		m_degree_limit = r;
		m_psyms = names2psyms(vs);
		m_mode = p_deg;
	}

	void degree::print(std::ostream &stream)
	{
		switch (m_mode) {
			case inactive:
				stream << "No degree limit set.";
				break;
			case deg:
				stream << "Degree limit: " << m_degree_limit;
				break;
			case p_deg:
				stream << "Partial degree limit: " << m_degree_limit << ", Affected symbols: [";
				for (std::size_t i = 0; i < m_psyms.size(); ++i) {
					stream << '\'' << m_psyms[i].get_name() << '\'';
					if (i < m_psyms.size() - 1) {
						stream << " ";
					}
				}
				stream << ']';
		}
	}
}
