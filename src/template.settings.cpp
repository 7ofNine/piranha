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

#include <boost/integer_traits.hpp>

#include "core/memory.h"
#include "core/settings.h"
#include "core/stream_manager.h"
#include "core/version.h"

namespace piranha
{
// Settings manager's static members.
	float settings::hash_max_load_factor = (float)(0.3);
	double settings::m_numerical_zero = 1E-80;
	const max_fast_uint settings::min_u = boost::integer_traits<max_fast_uint>::min();
	const max_fast_uint settings::max_u = boost::integer_traits<max_fast_uint>::max();
	const max_fast_int settings::min_i = boost::integer_traits<max_fast_int>::min();
	const max_fast_int settings::max_i = boost::integer_traits<max_fast_int>::max();
	std::string settings::m_path = "@THEORIES_INSTALL_PATH@";
	const std::string settings::m_default_path = "@THEORIES_INSTALL_PATH@";
	bool settings::m_debug = false;
	const std::string settings::m_version = __PIRANHA_VERSION;
	bool settings::enable_progress_display = true;
	settings::startup_class settings::startup;
#ifdef _PIRANHA_TBB
	const tbb::task_scheduler_init settings::tbb_init;
#endif

	settings::startup_class::startup_class()
	{
		// Startup report.
		std::cout << "This is Piranha version " << m_version << __PIRANHA_VERSION_TAG << std::endl;
		std::cout << "Default parameters initialized:" << std::endl;
		std::cout << "Print precision\t\t\t=\t" << stream_manager::digits() << std::endl;
		std::cout << "Numerical zero\t\t\t=\t" << numerical_zero() << std::endl;
		std::cout << "Fast unsigned int range\t\t=\t" << "[0," << max_u << ']' << std::endl;
		std::cout << "Fast signed int range\t\t=\t" << '[' << min_i << ',' << max_i << ']' << std::endl;
		std::cout << "Default path\t=\t" << m_default_path << std::endl;
		std::cout << "Piranha is ready." << std::endl;
		std::cout << "_______________________________" << std::endl << std::endl;
		// Setup cout.
		stream_manager::setup_print(std::cout);
#ifdef _PIRANHA_TBB
		std::cout << "TBB version strings:" << std::endl;
		std::cout << __TBB_VERSION_STRINGS;
		std::cout << "TBB:  BUILD_DATE\t\t" << __TBB_DATETIME << std::endl;
#endif
	}

/// Set path to theories of motion.
	void settings::set_path(const std::string &str)
	{
		m_path = str;
	}

/// Get version.
	const std::string &settings::version()
	{
		return m_version;
	}
}
