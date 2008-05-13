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

#include "core/stream_manager.h"

namespace piranha
{
	unsigned int stream_manager::m_digits = 15;
	stream_manager::out_format stream_manager::m_format = stream_manager::plain;
	const unsigned int stream_manager::m_min_digits = 0;
	const unsigned int stream_manager::m_max_digits = 50;
	stream_manager::fp_representation stream_manager::m_fp_rep = stream_manager::scientific;

	unsigned int stream_manager::digits()
	{
		return m_digits;
	}

	unsigned int stream_manager::min_digits()
	{
		return m_min_digits;
	}

	unsigned int stream_manager::max_digits()
	{
		return m_max_digits;
	}

	void stream_manager::set_digits(int n)
	{
		if (n < (int)m_min_digits || n > (int)m_max_digits) {
			std::cout << "Invalid number of digits." << std::endl;
		} else {
			m_digits = (unsigned int)n;
		}
	}

	void stream_manager::set_fp_rep(fp_representation fpr)
	{
		m_fp_rep = fpr;
	}

	stream_manager::fp_representation stream_manager::fp_rep()
	{
		return m_fp_rep;
	}

	void stream_manager::setup_print(std::ostream &out_stream)
	{
		out_stream << std::setprecision(m_digits);
		switch (m_fp_rep) {
		case scientific:
			out_stream << std::scientific;
			break;
		case decimal:
			out_stream << std::fixed;
			break;
		}
	}

	stream_manager::out_format stream_manager::format()
	{
		return m_format;
	}

	void stream_manager::set_format(out_format fmt)
	{
		m_format = fmt;
	}
}
