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



#include <cstddef>
#include <string>

#include "base_classes/base_counting_allocator.h"
#include "config.h"
#include "exceptions.h"
#include "math.h" // For static base-2 logarithm.

#pragma warning (push)
#pragma warning (disable: 4251)
namespace piranha
{
	/// Manager class for piranha-specific settings.
	/**
	 * This class manages parameters specific to piranha classes.
	 */
	class PIRANHA_VISIBLE settings
	{
		public:
			enum MultiplicationAlgorithm
			{
				ALGORITHM_AUTOMATIC    = 0,
				ALGORITHM_PLAIN        = 1,
				ALGORITHM_VECTOR_CODED = 2,
				ALGORITHM_HASH_CODED   = 3
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
				PIRANHA_THROW(not_implemented_error,"debug support was not compiled in");
#endif
			}


			/// Get Piranha version.
			static const char *get_version(); // no stl interface over dll border
			/// Cache size in kilobytes.
			static const std::size_t cache_size;// = _PIRANHA_CACHE_SIZE;
//			static_assert(cache_size > 0 && lg<cache_size>::value > 1, "Invalid value for cache size.");
			static bool blocker;
			static std::size_t get_max_pretty_print_size();
			static void set_max_pretty_print_size(int);
			static std::size_t get_nthread();
			static void set_nthread(const int &);

			static MultiplicationAlgorithm getMultiplicationAlgorithm()
			{
				return multiplicationAlgorithm;
			}


			static void setMultiplicationAlgorithm(MultiplicationAlgorithm algorithm)
			{
				if (algorithm < 0 || algorithm > 3)
                {
					PIRANHA_THROW(value_error, "Invalid multiplication algorithm");
				}
				multiplicationAlgorithm = algorithm;
			}


		private:

			// Startup class.
			class PIRANHA_VISIBLE startup_class
			{
				public:
					startup_class();
			};


			static std::size_t		        m_max_pretty_print_size;
			static std::size_t		        m_memory_limit;          // Memory limit in bytes.
			static double			        m_numerical_zero;        // Numerical zero.
			static bool			            m_debug;
			static const std::string	    m_version;
			static startup_class		    startup;
			static unsigned			        m_nthread;               // Number of threads available.
			static MultiplicationAlgorithm	multiplicationAlgorithm;
	};
}

#pragma warning (pop)
// Debug mode.
#ifndef NDEBUG
#define PIRANHA_DEBUG(statement) {if (settings::get_debug()) {statement;}}
#else
#define PIRANHA_DEBUG(statement)
#endif

#endif
