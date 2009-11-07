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

#ifndef PIRANHA_SETTINGS_H
#define PIRANHA_SETTINGS_H

#include <boost/integer_traits.hpp>
#include <cstddef>
#include <iostream>
#include <string>

#include "base_classes/base_counting_allocator.h"
#include "config.h"
#include "exceptions.h"
#include "integer_typedefs.h"
#include "math.h" // For static base-2 logarithm.

namespace piranha
{
	/// Manager class for piranha-specific settings.
	/**
	 * This class manages parameters specific to piranha classes.
	 */
	class __PIRANHA_VISIBLE settings
	{
		public:
			enum fp_representation {
				scientific,
				decimal
			};
			static std::size_t get_used_memory()
			{
				return base_counting_allocator::count();
			}
			static std::size_t get_memory_limit()
			{
				return m_memory_limit;
			}
			static void set_memory_limit(const std::size_t &limit)
			{
				m_memory_limit = limit;
			}
			// Getters.
			/// Get path to theories of motion files.
			static const std::string &get_path()
			{
				return m_path;
			}
			/// Get debug flag.
			static bool get_debug()
			{
				return m_debug;
			}
			static const double &get_numerical_zero()
			{
				return m_numerical_zero;
			}
			/// Set debug flag.
			static void set_debug(const bool &flag)
			{
#ifndef NDEBUG
				m_debug = flag;
#else
				(void)flag;
				piranha_throw(not_implemented_error,"debug support was not compiled in");
#endif
			}
			/// Get Piranha version.
			static const std::string &get_version();
			static void set_path(const std::string &);
			// Getters and setters for I/O.
			static std::size_t get_digits();
			static void set_digits(const int &n) {
				if (n < static_cast<int>(m_min_digits) || n > static_cast<int>(m_max_digits)) {
					piranha_throw(value_error,"invalid number of digits");
				} else {
					m_digits = static_cast<std::size_t>(n);
				}
			}
			static std::size_t get_min_digits();
			static std::size_t get_max_digits();
			static void setup_stream(std::ostream &);
			static fp_representation get_fp_repr();
			static void set_fp_repr(fp_representation);
			/// Cache size in kilobytes.
			static const std::size_t cache_size = _PIRANHA_CACHE_SIZE;
			p_static_check(cache_size > 0 && lg<cache_size>::value > 1, "Invalid value for cache size.");
			static bool blocker;
			static std::size_t get_max_pretty_print_size();
			static void set_max_pretty_print_size(int);
			static const std::size_t &get_nthread();
			static void set_nthread(const int &);
		private:
			/// Startup class.
			/**
			 * Startup class is constructed at piranha invocation and sets default parameters.
			 */
			class __PIRANHA_VISIBLE startup_class
			{
				public:
					startup_class();
			};
			static std::size_t		m_max_pretty_print_size;
			/// Memory limit in bytes.
			static std::size_t		m_memory_limit;
			/// Numerical zero.
			static double			m_numerical_zero;
			/// Path to theories of motion.
			static std::string		m_path;
			static std::string		m_default_path;
			static bool			m_debug;
			static const std::string	m_version;
			static startup_class		startup;
			/// Minimum number of digits for output streams.
			static const std::size_t	m_min_digits = 0;
			/// Maximum number of digits for output streams.
			static const std::size_t	m_max_digits = 50;
			/// Number of digits to display in output stream.
			static std::size_t		m_digits;
			/// Floating point representation.
			static fp_representation	m_fp_repr;
			/// Number of threads available.
			static size_t			m_nthread;
	};
}

// Debug mode.
#ifndef NDEBUG
#define __PDEBUG(statement) {if (settings::debug()) {statement;}}
#else
#define __PDEBUG(statement)
#endif

#endif
