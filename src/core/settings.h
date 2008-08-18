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

#include <iostream>
#include <string>

#include "base_classes/base_counting_allocator.h"
#include "config.h"
#include "exceptions.h"
#include "integer_typedefs.h"
#include "piranha_tbb.h" // For task scheduler init.

namespace piranha
{
	/// Manager class for piranha-specific settings.
	/**
	 * This class manages parameters specific to piranha classes.
	 */
	class __PIRANHA_VISIBLE settings
	{
		public:
			static size_t used_memory() {
				return base_counting_allocator::count();
			}
			static size_t memory_limit() {
				return m_memory_limit;
			}
			static void memory_limit(const size_t &limit) {
				m_memory_limit = limit;
			}
			/// Return maximum load factor for hashed containers.
			static const double &load_factor() {
				return m_hash_max_load_factor;
			}
			// Getters.
			/// Get numerical zero.
			static const double &numerical_zero() {
				return m_numerical_zero;
			}
			/// Get path to theories of motion files.
			static const std::string &path() {
				return m_path;
			}
			/// Get debug flag.
			static bool debug() {
				return m_debug;
			}
			/// Set debug flag.
			static void debug(const bool &flag) {
#ifdef _PIRANHA_DEBUG
				m_debug = flag;
#else
				(void)flag;
				throw(unsuitable("Debug support was not compiled in."));
#endif
			}
			/// Get Piranha version.
			static const std::string &version();
			/// Set maximum load factor for hashed containers.
			/**
			 * @see settings::load_factor().
			 */
			static void load_factor(const double &value) {
				if (value <= 0 || value >= 1) {
					throw(unsuitable("Please insert a real number in the ]0,1[ interval."));
				}
				m_hash_max_load_factor = value;
			}
			static void set_path(const std::string &);
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
			/// Memory limit in bytes.
			static size_t							m_memory_limit;
			/// Load factor for hashed containers.
			static double							m_hash_max_load_factor;
			/// Numerical zero.
			static double							m_numerical_zero;
			/// Minimum fast unsigned integer.
			static const max_fast_uint				min_u;
			/// Maximum fast unsigned integer.
			static const max_fast_uint				max_u;
			/// Minimum fast integer.
			static const max_fast_int				min_i;
			/// Maximum fast integer.
			static const max_fast_int				max_i;
			/// Path to theories of motion.
			static std::string						m_path;
			static const std::string				m_default_path;
			static bool								m_debug;
			static const std::string				m_version;
			static bool								enable_progress_display;
			static startup_class					startup;
#ifdef _PIRANHA_MT
			static const tbb::task_scheduler_init	tbb_init;
#endif
	};
}

// Debug mode.
#ifdef _PIRANHA_DEBUG
#define __PDEBUG(statement) {if (settings::debug()) {statement;}}
#else
#define __PDEBUG(statement)
#endif

#endif
