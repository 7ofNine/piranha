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
			enum out_format {
				plain,
				pretty,
				latex
			};
			enum fp_representation {
				scientific,
				decimal
			};
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
			// Getters and setters for I/O.
			static size_t digits();
			static void digits(const max_fast_int &n) {
				if (n < static_cast<max_fast_int>(m_min_digits) || n > static_cast<max_fast_int>(m_max_digits)) {
					throw unsuitable("Invalid number of digits.");
				} else {
					m_digits = static_cast<size_t>(n);
				}
			}
			static size_t min_digits();
			static size_t max_digits();
			static void setup_stream(std::ostream &);
			static out_format format();
			static void format(out_format);
			static fp_representation fp_repr();
			static void fp_repr(fp_representation);
			static void pi_simplify(const bool &);
			static bool pi_simplify();
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
			static const max_fast_uint				min_u = boost::integer_traits<max_fast_uint>::const_min;
			/// Maximum fast unsigned integer.
			static const max_fast_uint				max_u = boost::integer_traits<max_fast_uint>::const_max;
			/// Minimum fast integer.
			static const max_fast_int				min_i = boost::integer_traits<max_fast_int>::const_min;
			/// Maximum fast integer.
			static const max_fast_int				max_i = boost::integer_traits<max_fast_int>::const_max;
			/// Path to theories of motion.
			static std::string						m_path;
			static std::string						m_default_path;
			static bool								m_debug;
			static const std::string				m_version;
			static bool								enable_progress_display;
			static startup_class					startup;
			/// Minimum number of digits for output streams.
			static const size_t						m_min_digits = 0;
			/// Maximum number of digits for output streams.
			static const size_t						m_max_digits = 50;
			/// Number of digits to display in output stream.
			static size_t							m_digits;
			/// Format for output.
			static out_format						m_format;
			/// Floating point representation.
			static fp_representation				m_fp_repr;
			/// Simplify pi symbol in Fourier and Poisson series.
			static bool								m_pi_simplify;
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
