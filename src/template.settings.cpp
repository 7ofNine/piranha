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
#include <cstring>
#include <gmp.h>

#include "core/config.h"
#include "core/exceptions.h"
#include "core/integer_typedefs.h"
#include "core/memory.h"
#include "core/settings.h"
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
		memcpy(retval, static_cast<void const *>(ptr), std::min<size_t>(old_size,new_size));
		a.deallocate(static_cast<char *>(ptr),old_size);
		return retval;
	}

	static inline void gmp_free_func(void *ptr, size_t size) {
		std_counting_allocator<char> a;
		a.deallocate(static_cast<char *>(ptr),size);
	}

	// Settings' static members.
	size_t settings::m_memory_limit = 1500000000u; // ~ 1.5GByte
	double settings::m_hash_max_load_factor = 0.5;
	double settings::m_numerical_zero = 1E-80;
	// TODO: this one must be initialised to a better value in windows.
	const std::string settings::m_default_path = "@PIRANHA_INSTALL_PREFIX@/@THEORIES_INSTALL_PATH@";
	std::string settings::m_path = settings::m_default_path;
	bool settings::m_debug = false;
	const std::string settings::m_version = "@PIRANHA_VERSION@";
	bool settings::enable_progress_display = true;
	settings::startup_class settings::startup;
	size_t settings::m_digits = 15;
	settings::out_format settings::m_format = settings::plain;
	settings::fp_representation settings::m_fp_repr = settings::scientific;
#ifdef _PIRANHA_MT
	const tbb::task_scheduler_init settings::tbb_init;
#endif

	settings::startup_class::startup_class()
	{
		p_static_check(sizeof(char) == 1, "Wrong char size.");
		p_static_check(sizeof(char) == sizeof(bool), "Wrong char-bool size ratio.");
		// Startup report.
		std::cout << "Piranha version: " << m_version << '\n';
		std::cout << "Revision number: " << "@PIRANHA_REV_NUMBER@\n";
		std::cout << "Default parameters initialized:\n";
		std::cout << "Print precision = " << settings::digits() << '\n';
		std::cout << "Numerical zero = " << numerical_zero() << '\n';
		std::cout << "Fast unsigned int range = " << "[0," << max_u << ']' << '\n';
		std::cout << "Fast signed int range = " << '[' << min_i << ',' << max_i << ']' << '\n';
		std::cout << "Default path = " << m_default_path << '\n';
		std::cout << "Piranha is ready.\n";
		std::cout << "_______________________________" << '\n' << '\n';
		// Setup cout.
		settings::setup_stream(std::cout);
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

	size_t settings::digits()
	{
		return m_digits;
	}

	size_t settings::min_digits()
	{
		return m_min_digits;
	}

	size_t settings::max_digits()
	{
		return m_max_digits;
	}

	void settings::fp_repr(fp_representation fpr)
	{
		m_fp_repr = fpr;
	}

	settings::fp_representation settings::fp_repr()
	{
		return m_fp_repr;
	}

	void settings::setup_stream(std::ostream &out_stream)
	{
		out_stream << std::setprecision(m_digits);
		switch (m_fp_repr) {
		case scientific:
			out_stream << std::scientific;
			break;
		case decimal:
			out_stream << std::fixed;
			break;
		}
	}

	settings::out_format settings::format()
	{
		return m_format;
	}

	void settings::format(out_format fmt)
	{
		m_format = fmt;
	}
}
