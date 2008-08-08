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

#include <algorithm>
#include <boost/integer_traits.hpp>
#include <cstring>
#include <gmp.h>

#include "core/config.h"
#include "core/integer_typedefs.h"
#include "core/memory.h"
#include "core/settings.h"
#include "core/stream_manager.h"
#include "core/piranha_tbb.h" // For task scheduler init.

namespace piranha
{
	static inline void *gmp_alloc_func(size_t size) {
		std_counting_allocator<char> a;
		return static_cast<void *>(a.allocate(size));
	}

	static inline void *gmp_realloc_func(void *ptr, size_t old_size, size_t new_size) {
		std_counting_allocator<char> a;
		void *retval = static_cast<void *>(a.allocate(new_size));
		memcpy(retval, static_cast<void const *>(ptr), std::min(old_size,new_size));
		a.deallocate(static_cast<char *>(ptr),old_size);
		return retval;
	}

	static inline void gmp_free_func(void *ptr, size_t size) {
		std_counting_allocator<char> a;
		a.deallocate(static_cast<char *>(ptr),size);
	}

	// Settings manager's static members.
	max_fast_int settings::m_memory_limit = 1000000000; // ~ 1GByte
	float settings::hash_max_load_factor = (float)(0.3);
	double settings::m_numerical_zero = 1E-80;
	const max_fast_uint settings::min_u = boost::integer_traits<max_fast_uint>::min();
	const max_fast_uint settings::max_u = boost::integer_traits<max_fast_uint>::max();
	const max_fast_int settings::min_i = boost::integer_traits<max_fast_int>::min();
	const max_fast_int settings::max_i = boost::integer_traits<max_fast_int>::max();
	const std::string settings::m_default_path = "@PIRANHA_INSTALL_PREFIX@/@THEORIES_INSTALL_PATH@";
	std::string settings::m_path = settings::m_default_path;
	bool settings::m_debug = false;
	const std::string settings::m_version = "@PIRANHA_VERSION@";
	bool settings::enable_progress_display = true;
	settings::startup_class settings::startup;
#ifdef _PIRANHA_MT
	const tbb::task_scheduler_init settings::tbb_init;
#endif

	settings::startup_class::startup_class()
	{
		p_static_check(sizeof(char) == 1, "Wrong char size.");
		// Startup report.
		std::cout << "Piranha version: " << m_version << '\n';
		std::cout << "Revision number: " << "@PIRANHA_REV_NUMBER@\n";
		std::cout << "Default parameters initialized:\n";
		std::cout << "Print precision = " << stream_manager::digits() << '\n';
		std::cout << "Numerical zero = " << numerical_zero() << '\n';
		std::cout << "Fast unsigned int range = " << "[0," << max_u << ']' << '\n';
		std::cout << "Fast signed int range = " << '[' << min_i << ',' << max_i << ']' << '\n';
		std::cout << "Default path = " << m_default_path << '\n';
		std::cout << "Piranha is ready.\n";
		std::cout << "_______________________________" << '\n' << '\n';
		// Setup cout.
		stream_manager::setup_print(std::cout);
		// Setup GMP's memory allocation functions.
		mp_set_memory_functions(gmp_alloc_func, gmp_realloc_func, gmp_free_func);
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
