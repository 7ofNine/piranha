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
#include <boost/thread/thread.hpp>
#include <climits>
#include <cstddef>
#include <cstdlib> // For getenv in Windows.
#include <cstring>

#include "core/config.h"
#include "core/exceptions.h"
#include "core/integer_typedefs.h"
#include "core/memory.h"
#include "core/settings.h"

namespace piranha
{
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
	std::size_t settings::m_memory_limit = 1500000000u; // ~ 1.5GByte
	double settings::m_numerical_zero;
	std::string settings::m_default_path;
	std::string settings::m_path;
	bool settings::m_debug = false;
	const std::string settings::m_version = "@PIRANHA_VERSION@";
	std::size_t settings::m_digits;
	settings::fp_representation settings::m_fp_repr = settings::scientific;
	const std::size_t settings::cache_size;
	bool settings::blocker = false;
	settings::startup_class settings::startup;
	std::size_t settings::m_max_pretty_print_size = 500;
	std::size_t settings::m_nthread = 1;

	settings::startup_class::startup_class()
	{
		// Some static checks.
		// This check is here because in integer_typedefs we select integers based on the number of bits.
		p_static_check(CHAR_BIT == 8, "Invalid size of char.");
		// This is because somewhere we use chars as bools.
		p_static_check(sizeof(char) == sizeof(bool), "Wrong char-bool size ratio.");
		// TODO: where is this used? Find out and wipe it.
		p_static_check(sizeof(std::size_t) == sizeof(void *), "std::size_t and void * are not the same size.");
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
		// Setup number of threads.
		const std::size_t nthread = boost::thread::hardware_concurrency();
		if (!nthread) {
			std::cout << "Unable to detect automatically the number of hardware threads, setting value to 1.\n";
			set_nthread(1);
		} else {
			set_nthread(nthread);
		}
		// Startup report.
		std::cout << "Piranha version: " << m_version << '\n';
		std::cout << "Number of hardware threads: " << get_nthread() << '\n';
		std::cout << "Piranha GIT revision: " << "@PIRANHA_GIT_REVISION@" << '\n';
		std::cout << "Piranha is ready.\n";
		std::cout << "_______________________________" << '\n' << '\n';
		// Setup cout.
		settings::setup_stream(std::cout);
	}

	/// Set path to theories of motion.
	void settings::set_path(const std::string &str)
	{
		m_path = str;
	}

	/// Get version.
	const std::string &settings::get_version()
	{
		return m_version;
	}

	std::size_t settings::get_digits()
	{
		return m_digits;
	}

	std::size_t settings::get_min_digits()
	{
		return m_min_digits;
	}

	std::size_t settings::get_max_digits()
	{
		return m_max_digits;
	}

	void settings::set_fp_repr(fp_representation fpr)
	{
		m_fp_repr = fpr;
	}

	settings::fp_representation settings::get_fp_repr()
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

	std::size_t settings::get_max_pretty_print_size()
	{
		return m_max_pretty_print_size;
	}

	void settings::set_max_pretty_print_size(int n)
	{
		if (n < 10) {
			piranha_throw(value_error,"invalid max size for pretty printing, "
				"please insert an integer greater than 10");
		}
		m_max_pretty_print_size = n;
	}

	void settings::set_nthread(const int &n)
	{
		if (n <= 0) {
			piranha_throw(value_error,"invalid number of threads, "
				"please insert an integer greater than 0");
		}
		m_nthread = n;
	}

	const size_t &settings::get_nthread()
	{
		return m_nthread;
	}
}
