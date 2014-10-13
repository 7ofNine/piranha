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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/thread.hpp>
#include <climits>
#include <cstddef>

#include "core/config.h"
#include "core/exceptions.h"
#include "core/memory.h"
#include "core/settings.h"

namespace piranha
{
	// Settings' static members.
	std::size_t settings::m_memory_limit = 1500000000u; // ~ 1.5GByte
	double settings::m_numerical_zero = 1E-80;
	bool settings::m_debug = true;
	const std::string settings::m_version = "@PIRANHA_VERSION_STRING@";
	const std::size_t settings::cache_size = _PIRANHA_CACHE_SIZE;
	bool settings::blocker = false;
	unsigned settings::m_nthread = boost::thread::hardware_concurrency() ? boost::thread::hardware_concurrency() : 1;
	settings::startup_class settings::startup;
	std::size_t settings::m_max_pretty_print_size = 500;
	settings::multiplication_algorithm settings::m_mult_algo = settings::automatic;

	settings::startup_class::startup_class()
	{
		// Some static checks.
		// This check is here because in integer_typedefs we select integers based on the number of bits.
		PIRANHA_STATIC_CHECK(CHAR_BIT == 8, "Invalid size of char.");
		// This is because somewhere we use chars as bools.
		PIRANHA_STATIC_CHECK(sizeof(char) == sizeof(bool), "Wrong char-bool size ratio.");
		// TODO: where is this used? Find out and wipe it.
		PIRANHA_STATIC_CHECK(sizeof(std::size_t) == sizeof(void *), "std::size_t and void * are not the same size.");
		PIRANHA_STATIC_CHECK(__PIRANHA_MAX_ECHELON_LEVEL >= 0, "Max echelon level must be nonnegative.");
		PIRANHA_STATIC_CHECK(settings::cache_size > 0 && lg<settings::cache_size>::value > 1, "Invalid value for cache size.");
		// Startup report.
		std::cout << "Piranha version: " << "@PIRANHA_VERSION_STRING@" << '\n';
		std::cout << "Piranha GIT revision: " << "@PIRANHA_GIT_REVISION@" << '\n';
		std::cout << "Number of cores detected: " <<
			(boost::thread::hardware_concurrency() ? boost::thread::hardware_concurrency() : 1) << '\n';
		std::cout << "Piranha is ready.\n";
		std::cout << "_______________________________" << '\n' << '\n';
	}

	/// Get version.
	const std::string &settings::get_version()
	{
		return m_version;
	}

	std::size_t settings::get_max_pretty_print_size()
	{
		return m_max_pretty_print_size;
	}

	void settings::set_max_pretty_print_size(int n)
	{
		if (n < 10) {
			PIRANHA_THROW(value_error,"invalid max size for pretty printing, "
				"please insert an integer greater than 10");
		}
		m_max_pretty_print_size = boost::numeric_cast<std::size_t>(n);
	}

	void settings::set_nthread(const int &n)
	{
		if (n <= 0) {
			PIRANHA_THROW(value_error,"invalid number of threads, "
				"please insert an integer greater than 0");
		}
		m_nthread = boost::numeric_cast<unsigned>(n);
	}

	unsigned settings::get_nthread()
	{
		return m_nthread;
	}
}
