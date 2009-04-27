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
#include <boost/algorithm/string/replace.hpp> // For replacing "\" with "/" in Windows.
#include <cstring>
#include <cstdlib> // For getenv in Windows.
#include <gmp.h>

#include "core/config.h"
#include "core/exceptions.h"
#include "core/integer_typedefs.h"
#include "core/memory.h"
#include "core/settings.h"

namespace piranha
{
	extern "C" {
	static inline void *gmp_alloc_func(size_t size)
	{
		std_counting_allocator<char> a;
		return static_cast<void *>(a.allocate(size));
	}

	static inline void *gmp_realloc_func(void *ptr, size_t old_size, size_t new_size)
	{
		std_counting_allocator<char> a;
		void *retval = static_cast<void *>(a.allocate(new_size));
		memcpy(retval, static_cast<void const *>(ptr), std::min<size_t>(old_size,new_size));
		a.deallocate(static_cast<char *>(ptr),old_size);
		return retval;
	}

	static inline void gmp_free_func(void *ptr, size_t size)
	{
		std_counting_allocator<char> a;
		a.deallocate(static_cast<char *>(ptr),size);
	}
	}

	static inline std::string get_env_variable(const char *str)
	{
		if (getenv(str) != 0) {
			return std::string(getenv(str));
		} else {
			std::cout << "The environment variable '" << str << "' is not set.\n";
			return std::string();
		}
	}

	// Settings' static members.
	size_t settings::m_memory_limit = 1500000000u; // ~ 1.5GByte
	double settings::m_hash_max_load_factor = 1;
	double settings::m_numerical_zero;
	std::string settings::m_default_path;
	std::string settings::m_path;
	bool settings::m_debug = false;
	const std::string settings::m_version = "@PIRANHA_VERSION@";
	bool settings::enable_progress_display = true;
	size_t settings::m_digits;
	settings::out_format settings::m_format = settings::plain;
	settings::fp_representation settings::m_fp_repr = settings::scientific;
	const size_t settings::cache_size;
	bool settings::blocker = false;
	settings::startup_class settings::startup;
	size_t settings::m_max_pretty_print_size = 500;

	settings::startup_class::startup_class()
	{
		// Some static checks.
		p_static_check(sizeof(char) == 1, "Wrong char size.");
		p_static_check(sizeof(char) == sizeof(bool), "Wrong char-bool size ratio.");
		p_static_check(sizeof(max_fast_int) == sizeof(void *), "max_fast_int and void * are not the same size.");
		p_static_check(sizeof(size_t) == sizeof(void *), "size_t and void * are not the same size.");
		p_static_check(__PIRANHA_MAX_ECHELON_LEVEL >= 0, "Max echelon level must be nonnegative.");
		// Init values.
		m_numerical_zero = 1E-80;
		m_digits = 15;
		m_default_path =
#ifdef _PIRANHA_WIN32
			// NOTE: this is a bit hackish, but it works. The alternative would be to use
			// boost::filesystem but it looks like overkill, since it would be used just here.
			boost::algorithm::replace_all_copy(
				get_env_variable("ProgramFiles")+std::string("/Piranha @PIRANHA_VERSION@/@THEORIES_INSTALL_PATH@"),
				std::string("\\"),
				std::string("/")
			);
#else
			"@PIRANHA_INSTALL_PREFIX@/@THEORIES_INSTALL_PATH@";
#endif
		m_path = m_default_path;
		// Startup report.
		std::cout << "Piranha version: " << m_version << '\n';
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

	size_t settings::get_max_pretty_print_size()
	{
		return m_max_pretty_print_size;
	}

	void settings::set_max_pretty_print_size(int n)
	{
		if (n < 10) {
			throw unsuitable("Invalid max size for pretty printing.");
		}
		m_max_pretty_print_size = n;
	}
}
